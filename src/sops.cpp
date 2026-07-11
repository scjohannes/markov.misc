#include <cpp11.hpp>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace cpp11;

[[cpp11::register]] doubles cpp_markov_propagate(
    doubles_matrix<> initial,
    list transitions,
    integers non_absorb,
    integers absorb) {
  const int n = initial.nrow();
  const int states = initial.ncol();
  const int times = transitions.size() + 1;
  writable::doubles out(n * times * states);
  std::fill(out.begin(), out.end(), 0.0);
  out.attr("dim") = writable::integers({n, times, states});

  auto index = [n, times](int i, int t, int k) {
    return i + n * t + n * times * k;
  };
  for (int k = 0; k < states; ++k) {
    for (int i = 0; i < n; ++i) {
      out[index(i, 0, k)] = initial(i, k);
    }
  }

  for (int t = 1; t < times; ++t) {
    doubles_matrix<> trans(transitions[t - 1]);
    for (int origin_pos = 0; origin_pos < non_absorb.size(); ++origin_pos) {
      const int origin = non_absorb[origin_pos] - 1;
      const int row_offset = origin_pos * n;
      for (int i = 0; i < n; ++i) {
        const double previous = out[index(i, t - 1, origin)];
        if (previous == 0.0) continue;
        for (int target = 0; target < states; ++target) {
          out[index(i, t, target)] += previous * trans(row_offset + i, target);
        }
      }
    }
    for (int a_pos = 0; a_pos < absorb.size(); ++a_pos) {
      const int a = absorb[a_pos] - 1;
      for (int i = 0; i < n; ++i) {
        out[index(i, t, a)] += out[index(i, t - 1, a)];
      }
    }
  }
  return out;
}

[[cpp11::register]] doubles cpp_normalize_probability_array(
    doubles values,
    int draws,
    int observations,
    int states) {
  writable::doubles out(values);
  auto index = [draws, observations](int d, int i, int k) {
    return d + draws * i + draws * observations * k;
  };
  for (int d = 0; d < draws; ++d) {
    for (int i = 0; i < observations; ++i) {
      double total = 0.0;
      bool finite = true;
      for (int k = 0; k < states; ++k) {
        const double value = out[index(d, i, k)];
        if (ISNAN(value)) continue;
        if (!std::isfinite(value)) finite = false;
        total += value;
      }
      if (finite && total > 0.0) {
        for (int k = 0; k < states; ++k) {
          const int pos = index(d, i, k);
          if (!ISNAN(out[pos])) out[pos] /= total;
        }
      }
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
  writable::doubles out(draws * group_count * times * states);
  std::fill(out.begin(), out.end(), 0.0);
  std::vector<int> counts(group_count, 0);
  auto input_index = [draws, observations, times](int d, int i, int t, int k) {
    return d + draws * i + draws * observations * t +
      draws * observations * times * k;
  };
  auto output_index = [draws, group_count, times](int d, int g, int t, int k) {
    return d + draws * g + draws * group_count * t +
      draws * group_count * times * k;
  };

  for (int i = 0; i < observations; ++i) {
    const int g = groups[i] - 1;
    if (g < 0 || g >= group_count) continue;
    counts[g] += 1;
    for (int k = 0; k < states; ++k) {
      for (int t = 0; t < times; ++t) {
        for (int d = 0; d < draws; ++d) {
          out[output_index(d, g, t, k)] += values[input_index(d, i, t, k)];
        }
      }
    }
  }
  for (int g = 0; g < group_count; ++g) {
    if (counts[g] == 0) continue;
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
  writable::doubles out(draws * observations * states);
  std::fill(out.begin(), out.end(), 0.0);
  auto state_index = [draws, observations](int d, int i, int k) {
    return d + draws * i + draws * observations * k;
  };
  const int origins = non_absorb.size();
  auto trans_index = [draws, observations, origins](
      int d, int i, int origin_pos, int target) {
    const int rows = observations * origins;
    return d + draws * (origin_pos * observations + i) + draws * rows * target;
  };

  for (int origin_pos = 0; origin_pos < non_absorb.size(); ++origin_pos) {
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
  }
  for (int a_pos = 0; a_pos < absorb.size(); ++a_pos) {
    const int a = absorb[a_pos] - 1;
    for (int i = 0; i < observations; ++i) {
      for (int d = 0; d < draws; ++d) {
        out[state_index(d, i, a)] += previous[state_index(d, i, a)];
      }
    }
  }
  out.attr("dim") = writable::integers({draws, observations, states});
  return out;
}
