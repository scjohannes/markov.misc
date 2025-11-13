#' Use bootstrap to compute confidence intervals for the time alive and out of 
#' hospital 
#'
#' Details...
#'
#' Output is the sample estimate of time alive and out of hospital, and an upper 
#' and lower confidence bound for the estimate.  
#'
#' @param data  A data frame containing the patient data
#' @param n_boot Number of bootstrap samples
#' @param formula Model formula for orm
#' @param workers Number of workers used for parallelization. Default is availableCores()-1
#' @param parallel Whether parallelization should be used
#' @param ylevels States in the data
#' @param absorb Absorbing state
#' @param times Time points in the data. Default is the maximum time observed in the data
#' @param varnames List of variable names in the data
#'
#' @keywords bootstrap time alive and out of hospital
#'
#' @importFrom rms orm
#' @importFrom rms soprobMarkovOrdm
#' @importFrom rsample group_bootstraps
#' @importFrom rsample analysis
#' @importFrom furrr plan
#' @importFrom furrr availableCores
#' @importFrom furrr future_map
#'
#' @examples
#' taooh_bootstraps()
#'
#' @export

# Planned extensions: 
# - If beta == TRUE, stop the function and return beta treatment (for beta bootstrapping)
# - Handle missing states in the data (change ylevels and absorb depanding on the data)
# - Use different bootstrapping procedure. 
#     - Up to now, some resamples contain less unique IDs than the original data. 
#       Unclear, why this can happen. 

taooh_bootstrap <- function(
    data, 
    n_boot, 
    formula,
    workers = NULL,
    parallel = TRUE,
    ylevels = 1:6,
    absorb = 6,
    times = 1:max(data[["time"]]),
    varnames = list(tvarname = "time", 
                     pvarname = "yprev", 
                     id = "id", 
                     tx = "tx")
    ) {
  
  # Function to compute time alive and out of hospital
  taooh <- function(data, formula, times, ylevels = 1:6, absorb = 6,
                    varnames = list(tvarname = "time", 
                                    pvarname = "yprev", 
                                    id = "id", 
                                    tx = "tx")
  ) {
    # Convert splits object to data.frame
    data <- analysis(data)
    
    # Define new id variable (unique to each bootstrap draw)
    data$block_count <- ave(data[[varnames$id]], 
                            data[[varnames$id]], data[[varnames$tvarname]], 
                            FUN = seq_along)
    data$new_id <- factor(paste(data[[varnames$id]], data$block_count, sep = "_"))
    
    # 1. Fit model
    m <- orm(formula, data = data)
    
    # 2. Generate covariate data.frame to predict on
    X <- data[!duplicated(data$new_id), ] # first row of each individual
    
    # 2. Run soprobMarkovOrdm to get state probability predictions of each individual
    sop_mat <- matrix(nrow = max(times), ncol = length(unique(data$new_id)))
    for (i in 1:length(unique(data$new_id))) {
      sop_mat[, i] <- soprobMarkovOrdm(
        object = m,
        data = X[i, ],
        times = times,
        ylevels = ylevels,
        absorb = absorb,
        tvarname = varnames$tvarname,
        pvarname = varnames$pvarname,
        gap = 1
      )[,1] # only save state 1 SOPs
    }
    colnames(sop_mat) <- unique(data$new_id)
    
    # 3. Indicate treatment and control ids: First col := IDs, Second col := Treatment indicator
    id_tx <- aggregate(data[[varnames$tx]], 
                       by = list(data$new_id), 
                       FUN = unique, simplify = TRUE)
    
    # 4. Compute TAOOH estimate by treatment group
    SOP_tx <- sum(rowMeans(sop_mat[, colnames(sop_mat) %in% id_tx[id_tx[,2] == 1,1]]))
    SOP_ctrl <- sum(rowMeans(sop_mat[, colnames(sop_mat) %in% id_tx[id_tx[,2] == 0,1]]))
    
    return(c(delta_taooh = SOP_tx - SOP_ctrl, SOP_tx = SOP_tx, SOP_ctrl = SOP_ctrl))
  }
  
  
  # Bootstrap samples
  resample <- group_bootstraps(data, 
                               group = id, times = n_boot, apparent = FALSE)
  
  # Arguments
  formula <- as.formula(formula)
  if(is.null(workers)) workers <- availableCores()-1
  
  # Apply the taooh() function to the bootstrap samples in parallel.
  if(parallel) {
    plan(multisession, workers = workers)
  } else {plan(sequential)}
  
  bs_SOP <- resample %>% 
    mutate(models = future_map(splits, 
                               \(.x) taooh(
                                 .x,
                                 formula = formula, 
                                 times = times, 
                                 ylevels = ylevels, 
                                 absorb = absorb,
                                 varnames = varnames)[1]
    ))
  
  plan(sequential)
  
  return(bs_SOP)
}
