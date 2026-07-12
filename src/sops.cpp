#include <cpp11.hpp>
#include <R_ext/Utils.h>
#include <R_ext/Random.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace cpp11;

namespace {

R_xlen_t checked_product(std::initializer_list<int> dimensions) {
  R_xlen_t result = 1;
  for (const int dimension : dimensions) {
    if (dimension < 0) {
      cpp11::stop("Native SOP dimensions must be non-negative.");
    }
    if (dimension != 0 &&
        result > R_XLEN_T_MAX / static_cast<R_xlen_t>(dimension)) {
      cpp11::stop("Native SOP dimensions exceed R's vector-size limit.");
    }
    result *= static_cast<R_xlen_t>(dimension);
  }
  return result;
}

void validate_state_partition(
    const integers& non_absorb,
    const integers& absorb,
    const int states) {
  std::vector<bool> seen(states, false);
  for (R_xlen_t pos = 0; pos < non_absorb.size(); ++pos) {
    const int state = non_absorb[pos];
    if (state == NA_INTEGER || state < 1 || state > states || seen[state - 1]) {
      cpp11::stop("Invalid or duplicate non-absorbing state index.");
    }
    seen[state - 1] = true;
  }
  for (R_xlen_t pos = 0; pos < absorb.size(); ++pos) {
    const int state = absorb[pos];
    if (state == NA_INTEGER || state < 1 || state > states || seen[state - 1]) {
      cpp11::stop("Invalid, duplicate, or overlapping absorbing state index.");
    }
    seen[state - 1] = true;
  }
  if (std::find(seen.begin(), seen.end(), false) != seen.end()) {
    cpp11::stop("Native SOP state indices must form a complete partition.");
  }
}

double stable_logistic(const double value) {
  if (ISNAN(value)) {
    return NA_REAL;
  }
  if (value >= 0.0) {
    return 1.0 / (1.0 + std::exp(-value));
  }
  const double exponent = std::exp(value);
  return exponent / (1.0 + exponent);
}

}  // namespace

[[cpp11::register]] doubles cpp_markov_propagate(
    doubles_matrix<> initial,
    list transitions,
    integers non_absorb,
    integers absorb) {
  const int n = initial.nrow();
  const int states = initial.ncol();
  const int times = transitions.size() + 1;
  validate_state_partition(non_absorb, absorb, states);
  writable::doubles out(checked_product({n, times, states}));
  std::fill(out.begin(), out.end(), 0.0);
  out.attr("dim") = writable::integers({n, times, states});

  auto index = [n, times](const int i, const int t, const int k) -> R_xlen_t {
    return static_cast<R_xlen_t>(i) +
      static_cast<R_xlen_t>(n) * t +
      static_cast<R_xlen_t>(n) * times * k;
  };
  for (int k = 0; k < states; ++k) {
    for (int i = 0; i < n; ++i) {
      out[index(i, 0, k)] = initial(i, k);
    }
  }

  for (int t = 1; t < times; ++t) {
    doubles_matrix<> trans(transitions[t - 1]);
    if (trans.nrow() != n * non_absorb.size() || trans.ncol() != states) {
      cpp11::stop("Transition matrix dimensions do not match the SOP state plan.");
    }
    for (R_xlen_t origin_pos = 0; origin_pos < non_absorb.size(); ++origin_pos) {
      const int origin = non_absorb[origin_pos] - 1;
      const R_xlen_t row_offset = origin_pos * n;
      for (int i = 0; i < n; ++i) {
        const double previous = out[index(i, t - 1, origin)];
        if (previous == 0.0) {
          continue;
        }
        for (int target = 0; target < states; ++target) {
          out[index(i, t, target)] +=
            previous * trans(row_offset + i, target);
        }
      }
    }
    for (R_xlen_t a_pos = 0; a_pos < absorb.size(); ++a_pos) {
      const int state = absorb[a_pos] - 1;
      for (int i = 0; i < n; ++i) {
        out[index(i, t, state)] += out[index(i, t - 1, state)];
      }
    }
    if ((t & 7) == 0) {
      cpp11::check_user_interrupt();
    }
  }
  return out;
}

[[cpp11::register]] doubles cpp_markov_update_logits(
    doubles_matrix<> previous,
    doubles_matrix<> logits,
    integers non_absorb,
    integers absorb) {
  const int n = previous.nrow();
  const int states = previous.ncol();
  const int thresholds = states - 1;
  const int origins = non_absorb.size();
  validate_state_partition(non_absorb, absorb, states);
  if (logits.nrow() != n * origins || logits.ncol() != thresholds) {
    cpp11::stop("Logit workspace dimensions do not match the SOP state plan.");
  }

  writable::doubles out(checked_product({n, states}));
  std::fill(out.begin(), out.end(), 0.0);
  out.attr("dim") = writable::integers({n, states});
  auto output_index = [n](const int i, const int state) -> R_xlen_t {
    return static_cast<R_xlen_t>(i) + static_cast<R_xlen_t>(n) * state;
  };
  std::vector<double> probabilities(states, 0.0);

  for (int origin_pos = 0; origin_pos < origins; ++origin_pos) {
    const int origin = non_absorb[origin_pos] - 1;
    const int row_offset = origin_pos * n;
    for (int i = 0; i < n; ++i) {
      const double mass = previous(i, origin);
      if (mass == 0.0) {
        continue;
      }
      const int row = row_offset + i;
      double prior_cumulative = stable_logistic(logits(row, 0));
      probabilities[0] = std::max(0.0, 1.0 - prior_cumulative);
      double total = probabilities[0];
      for (int threshold = 1; threshold < thresholds; ++threshold) {
        const double cumulative = stable_logistic(logits(row, threshold));
        probabilities[threshold] =
          std::max(0.0, prior_cumulative - cumulative);
        total += probabilities[threshold];
        prior_cumulative = cumulative;
      }
      probabilities[states - 1] = std::max(0.0, prior_cumulative);
      total += probabilities[states - 1];
      if (!std::isfinite(total) || total <= 0.0) {
        cpp11::stop("Invalid transition probabilities in native SOP propagation.");
      }
      const double scale = mass / total;
      for (int target = 0; target < states; ++target) {
        out[output_index(i, target)] += probabilities[target] * scale;
      }
    }
    if ((origin_pos & 7) == 0) {
      cpp11::check_user_interrupt();
    }
  }

  for (R_xlen_t a_pos = 0; a_pos < absorb.size(); ++a_pos) {
    const int state = absorb[a_pos] - 1;
    for (int i = 0; i < n; ++i) {
      out[output_index(i, state)] += previous(i, state);
    }
  }
  return out;
}

[[cpp11::register]] doubles cpp_markov_update_po(
    doubles_matrix<> previous,
    doubles scalar_predictor,
    doubles cutpoints,
    integers non_absorb,
    integers absorb) {
  const int n = previous.nrow();
  const int states = previous.ncol();
  const int thresholds = states - 1;
  const int origins = non_absorb.size();
  validate_state_partition(non_absorb, absorb, states);
  if (scalar_predictor.size() != checked_product({n, origins}) ||
      cutpoints.size() != thresholds) {
    cpp11::stop("PO predictor dimensions do not match the SOP state plan.");
  }

  writable::doubles out(checked_product({n, states}));
  std::fill(out.begin(), out.end(), 0.0);
  out.attr("dim") = writable::integers({n, states});
  auto output_index = [n](const int i, const int state) -> R_xlen_t {
    return static_cast<R_xlen_t>(i) + static_cast<R_xlen_t>(n) * state;
  };
  std::vector<double> probabilities(states, 0.0);

  for (int origin_pos = 0; origin_pos < origins; ++origin_pos) {
    const int origin = non_absorb[origin_pos] - 1;
    const R_xlen_t row_offset = static_cast<R_xlen_t>(origin_pos) * n;
    for (int i = 0; i < n; ++i) {
      const double mass = previous(i, origin);
      if (mass == 0.0) {
        continue;
      }
      const double scalar = scalar_predictor[row_offset + i];
      double prior_cumulative = stable_logistic(cutpoints[0] + scalar);
      probabilities[0] = std::max(0.0, 1.0 - prior_cumulative);
      double total = probabilities[0];
      for (int threshold = 1; threshold < thresholds; ++threshold) {
        const double cumulative = stable_logistic(cutpoints[threshold] + scalar);
        probabilities[threshold] =
          std::max(0.0, prior_cumulative - cumulative);
        total += probabilities[threshold];
        prior_cumulative = cumulative;
      }
      probabilities[states - 1] = std::max(0.0, prior_cumulative);
      total += probabilities[states - 1];
      if (!std::isfinite(total) || total <= 0.0) {
        cpp11::stop("Invalid transition probabilities in native PO propagation.");
      }
      const double scale = mass / total;
      for (int target = 0; target < states; ++target) {
        out[output_index(i, target)] += probabilities[target] * scale;
      }
    }
    if ((origin_pos & 7) == 0) {
      cpp11::check_user_interrupt();
    }
  }

  for (R_xlen_t a_pos = 0; a_pos < absorb.size(); ++a_pos) {
    const int state = absorb[a_pos] - 1;
    for (int i = 0; i < n; ++i) {
      out[output_index(i, state)] += previous(i, state);
    }
  }
  return out;
}

[[cpp11::register]] integers cpp_sample_categorical_rows(
    doubles_matrix<> probabilities) {
  const int observations = probabilities.nrow();
  const int categories = probabilities.ncol();
  if (categories < 1) {
    cpp11::stop("Categorical sampling requires at least one category.");
  }
  std::vector<double> totals(observations, 0.0);
  for (int i = 0; i < observations; ++i) {
    double total = 0.0;
    for (int category = 0; category < categories; ++category) {
      const double probability = probabilities(i, category);
      if (!std::isfinite(probability) || probability < 0.0) {
        cpp11::stop("Categorical probabilities must be finite and non-negative.");
      }
      total += probability;
    }
    if (!std::isfinite(total) || total <= 0.0) {
      cpp11::stop("Categorical probability rows must have positive mass.");
    }
    totals[i] = total;
  }

  writable::integers out(observations);
  GetRNGstate();
  for (int i = 0; i < observations; ++i) {
    const double draw = unif_rand() * totals[i];
    double cumulative = 0.0;
    int selected = categories;
    for (int category = 0; category < categories; ++category) {
      cumulative += probabilities(i, category);
      if (draw <= cumulative) {
        selected = category + 1;
        break;
      }
    }
    out[i] = selected;
  }
  PutRNGstate();
  return out;
}

[[cpp11::register]] doubles cpp_normalize_probability_array(
    doubles values,
    int draws,
    int observations,
    int states) {
  const R_xlen_t expected = checked_product({draws, observations, states});
  if (values.size() != expected) {
    cpp11::stop("Probability array length does not match its dimensions.");
  }
  writable::doubles out(values);
  auto index = [draws, observations](const int d, const int i, const int k) -> R_xlen_t {
    return static_cast<R_xlen_t>(d) +
      static_cast<R_xlen_t>(draws) * i +
      static_cast<R_xlen_t>(draws) * observations * k;
  };
  for (int d = 0; d < draws; ++d) {
    for (int i = 0; i < observations; ++i) {
      double total = 0.0;
      bool finite = true;
      for (int k = 0; k < states; ++k) {
        const double value = out[index(d, i, k)];
        if (ISNAN(value)) {
          continue;
        }
        if (!std::isfinite(value)) {
          finite = false;
        }
        total += value;
      }
      if (finite && total > 0.0) {
        for (int k = 0; k < states; ++k) {
          const R_xlen_t pos = index(d, i, k);
          if (!ISNAN(out[pos])) {
            out[pos] /= total;
          }
        }
      }
    }
    if ((d & 31) == 0) {
      cpp11::check_user_interrupt();
    }
  }
  out.attr("dim") = writable::integers({draws, observations, states});
  return out;
}

[[cpp11::register]] doubles cpp_reduce_sops_draws(
    doubles values,
    int draws,
    int observations,
    int times,
    int states,
    integers groups,
    int group_count) {
  const R_xlen_t expected = checked_product({draws, observations, times, states});
  if (values.size() != expected || groups.size() != observations || group_count < 1) {
    cpp11::stop("Posterior SOP dimensions and grouping indices are not aligned.");
  }
  writable::doubles out(checked_product({draws, group_count, times, states}));
  std::fill(out.begin(), out.end(), 0.0);
  std::vector<int> counts(group_count, 0);
  auto input_index = [draws, observations, times](
      const int d, const int i, const int t, const int k) -> R_xlen_t {
    return static_cast<R_xlen_t>(d) +
      static_cast<R_xlen_t>(draws) * i +
      static_cast<R_xlen_t>(draws) * observations * t +
      static_cast<R_xlen_t>(draws) * observations * times * k;
  };
  auto output_index = [draws, group_count, times](
      const int d, const int g, const int t, const int k) -> R_xlen_t {
    return static_cast<R_xlen_t>(d) +
      static_cast<R_xlen_t>(draws) * g +
      static_cast<R_xlen_t>(draws) * group_count * t +
      static_cast<R_xlen_t>(draws) * group_count * times * k;
  };

  for (int i = 0; i < observations; ++i) {
    const int group = groups[i];
    if (group == NA_INTEGER || group < 1 || group > group_count) {
      cpp11::stop("Invalid or missing SOP grouping index.");
    }
    const int g = group - 1;
    counts[g] += 1;
    for (int k = 0; k < states; ++k) {
      for (int t = 0; t < times; ++t) {
        for (int d = 0; d < draws; ++d) {
          out[output_index(d, g, t, k)] += values[input_index(d, i, t, k)];
        }
      }
    }
    if ((i & 255) == 0) {
      cpp11::check_user_interrupt();
    }
  }
  for (int g = 0; g < group_count; ++g) {
    if (counts[g] == 0) {
      continue;
    }
    const double scale = 1.0 / counts[g];
    for (int k = 0; k < states; ++k) {
      for (int t = 0; t < times; ++t) {
        for (int d = 0; d < draws; ++d) {
          out[output_index(d, g, t, k)] *= scale;
        }
      }
    }
  }
  out.attr("dim") = writable::integers({draws, group_count, times, states});
  return out;
}

[[cpp11::register]] doubles cpp_markov_update_draws(
    doubles previous,
    doubles transition,
    int draws,
    int observations,
    int states,
    integers non_absorb,
    integers absorb) {
  const R_xlen_t state_size = checked_product({draws, observations, states});
  const int origins = non_absorb.size();
  const R_xlen_t transition_size =
    checked_product({draws, observations, origins, states});
  if (previous.size() != state_size || transition.size() != transition_size) {
    cpp11::stop("Posterior transition dimensions do not match the SOP state plan.");
  }
  validate_state_partition(non_absorb, absorb, states);
  writable::doubles out(state_size);
  std::fill(out.begin(), out.end(), 0.0);
  auto state_index = [draws, observations](
      const int d, const int i, const int k) -> R_xlen_t {
    return static_cast<R_xlen_t>(d) +
      static_cast<R_xlen_t>(draws) * i +
      static_cast<R_xlen_t>(draws) * observations * k;
  };
  auto trans_index = [draws, observations, origins](
      const int d, const int i, const int origin_pos, const int target) -> R_xlen_t {
    const R_xlen_t rows = static_cast<R_xlen_t>(observations) * origins;
    return static_cast<R_xlen_t>(d) +
      static_cast<R_xlen_t>(draws) * (origin_pos * observations + i) +
      static_cast<R_xlen_t>(draws) * rows * target;
  };

  for (int origin_pos = 0; origin_pos < origins; ++origin_pos) {
    const int origin = non_absorb[origin_pos] - 1;
    for (int target = 0; target < states; ++target) {
      for (int i = 0; i < observations; ++i) {
        for (int d = 0; d < draws; ++d) {
          out[state_index(d, i, target)] +=
            previous[state_index(d, i, origin)] *
            transition[trans_index(d, i, origin_pos, target)];
        }
      }
    }
    if ((origin_pos & 7) == 0) {
      cpp11::check_user_interrupt();
    }
  }
  for (R_xlen_t a_pos = 0; a_pos < absorb.size(); ++a_pos) {
    const int state = absorb[a_pos] - 1;
    for (int i = 0; i < observations; ++i) {
      for (int d = 0; d < draws; ++d) {
        out[state_index(d, i, state)] += previous[state_index(d, i, state)];
      }
    }
  }
  out.attr("dim") = writable::integers({draws, observations, states});
  return out;
}
