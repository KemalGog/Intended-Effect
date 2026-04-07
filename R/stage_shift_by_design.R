##############################################################################
# INCIDENCE PROJECTION FUNCTIONS FOR STANDARD AND INTENDED EFFECT TRIAL DESIGNS
##############################################################################
# Authors: Jane Lange, Kemal Gogebakan, 3/19/2026
#
# Description:
# These functions implement a continuous-time Markov chain and HMM
# framework to project cancer detection in randomized cancer screening trials.
#
# Given trial design, natural history parameters, and test
# sensitivities, the functions estimate stage-specific cancer incidence in the
# control and screen arms over time.
#
# The framework supports both:
# - Standard design
# - Intended-effect (IE) design
#
# Outputs include:
# - Early-stage, late-stage, and overall incidence
# - Screen-detected and interval-detected cancers
# - Cancers detected in the control arm
# - Counts (scaled by arm size) and probabilities (incidence)
# - Stage shift (reduction in late-stage incidence)
#
# Results can be returned as yearly or cumulative measures, with optional
# inclusion of early-stage and overall outcomes.
#####################################################################
#####################################################################
# 1. NATURAL HISTORY MODEL FUNCTIONS
#####################################################################

#' Compute initial state distribution at starting age
#'
#' Computes the initial state distribution at the specified starting age,
#' conditional on not yet being clinically diagnosed. The distribution is
#' restricted to preclinical states (including healthy and preclinical disease).
#'
#' @param rate_matrix Transition rate matrix.
#' @param start_age Numeric value specifying the starting age.
#'
#' @return A numeric vector representing the initial state distribution over
#'   non-clinical states.
#' @export  
get_init <- function(rate_matrix, start_age) {
 
  k <- nrow(rate_matrix)
  init <- vector()
  init[k] <- 0
  init[k - 1] <- 0
  
  prob_mat_a1 <- MatrixExp(t = start_age, mat = rate_matrix)
  denom <- 1 - prob_mat_a1[1, (k - 1)] - prob_mat_a1[1, k]
  
  for (j in 1:(k - 2)) {
    init[j] <- prob_mat_a1[1, j] / denom
  }
  
  init
}

#' Construct emission probability matrix
#'
#' Builds the emission probability matrix linking latent disease states to
#' observed screen-detected and clinically detected states. Assumes perfect
#' specificity of tests.
#'
#' @param k Integer specifying the number of latent states.
#' @param sens_e Numeric value for early-stage sensitivity.
#' @param sens_l Numeric value for late-stage sensitivity.
#'
#' @return A matrix of emission probabilities mapping latent states to observed states.
get_emission <- function(k, sens_e, sens_l) {
  emit <- matrix(0, nrow = k, ncol = 5)
  
  emit[1:(k - 4), 1] <- 1
  
  emit[(k - 3), 2] <- sens_e
  emit[(k - 3), 1] <- 1 - sens_e
  
  emit[(k - 2), 3] <- sens_l
  emit[(k - 2), 1] <- 1 - sens_l
  
  emit[(k - 1), 4] <- 1
  emit[k, 5] <- 1
  
  emit
}

#' Create event time list for detection processes
#'
#' Creates a list of times corresponding to screen detection, interval detection,
#' or clinical detection in the post-screen period.
#'
#' @param screen_times Numeric vector of screening times.
#' @param preclin Logical indicating whether the event is screen-detected.
#' @param post Logical indicating whether the event is post-screen clinical detection.
#' @param extended_followup Numeric value specifying additional follow-up time after the last screen.
#'
#' @return A list of sequences of times for use in the HMM likelihood.
get_times <- function(screen_times, preclin = FALSE, post = FALSE, extended_followup = NULL) {

    if (post) {
    if (is.null(extended_followup)) stop("Need extended_followup when post = TRUE")
    return(list(c(screen_times, max(screen_times) + extended_followup)))
  }
  
  n <- length(screen_times)
  the_times <- list()
  
  if (n == 1) {
    if (preclin) {
      the_times[[1]] <- c(0)
    } else {
      if (is.null(extended_followup)) {
        stop("Need extended_followup after single exam when preclin = FALSE")
      } else {
        the_times[[1]] <- c(0, extended_followup)
      }
    }
    return(the_times)
  }
  
  if (preclin) {
    for (i in seq_len(n)) {
      the_times[[i]] <- screen_times[1:i]
    }
  } else {
    for (i in seq_len(n - 1)) {
      the_times[[i]] <- screen_times[1:(i + 1)]
    }
    if (!is.null(extended_followup)) {
      the_times[[length(the_times) + 1]] <- c(screen_times, tail(screen_times, 1) + extended_followup)
    }
  }
  
  the_times
}

#' Create observed state sequences for detection outcomes
#'
#' Creates observed state sequences corresponding to screen-detected,
#' interval-detected, or post-screen detected outcomes.
#'
#' @param screen_times Numeric vector of screening times.
#' @param stage Character string indicating stage ("early" or "late").
#' @param preclin Logical indicating whether the event is screen-detected.
#' @param post Logical indicating whether the event is post-screen clinical detection.
#' @param extended_followup Numeric value specifying additional follow-up time after the last screen.
#'
#' @return A list of sequences of observed states for use in the HMM computation.
get_data <- function(screen_times,
                     stage = c("early", "late"),
                     preclin = FALSE,
                     post = FALSE,
                     extended_followup = NULL) {

  stage <- match.arg(stage)
  
  preclin_code <- if (stage == "early") 2 else 3
  clinical_code <- if (stage == "early") 4 else 5
  
  if (post) {
    if (is.null(extended_followup)) stop("Need extended_followup when post = TRUE")
    return(list(c(rep(1, length(screen_times)), clinical_code)))
  }
  
  n <- length(screen_times)
  the_data <- list()
  
  if (n == 1) {
    if (preclin) {
      the_data[[1]] <- c(preclin_code)
    } else {
      if (is.null(extended_followup)) {
        stop("Need extended_followup after single exam when preclin = FALSE")
      } else {
        the_data[[1]] <- c(1, clinical_code)
      }
    }
    return(the_data)
  }
  
  if (preclin) {
    the_data[[1]] <- c(preclin_code)
    for (i in seq_len(n - 1)) {
      the_data[[i + 1]] <- c(rep(1, i), preclin_code)
    }
  } else {
    for (i in seq_len(n - 1)) {
      the_data[[i]] <- c(rep(1, i), clinical_code)
    }
    if (!is.null(extended_followup)) {
      the_data[[length(the_data) + 1]] <- c(rep(1, n), clinical_code)
    }
  }
  
  the_data
}

#' Compute transition probability matrices over observation intervals
#'
#' Computes lists of transition probability matrices for multiple observation
#' intervals under a time-homogeneous continuous-time Markov chain:
#' P(t2 - t1), P(t3 - t2), ..., P(tn - t_{n-1}).
#'
#' If an interval corresponds to an exact event time, that interval is evaluated
#' using a transition density rather than a transition probability. This is
#' implemented by multiplying the transition probability matrix by the
#' off-diagonal part of the rate matrix.
#'
#' @param time_intervals_list List of numeric vectors representing interval endpoints.
#' @param rate_matrix_list List of transition rate matrices.
#' @param exact_time_rank Optional integer index of the interval to evaluate as a density
#'   (i.e., exact event time).
#'
#' @return A list of lists of transition probability matrices, one for each observation pattern.
transition_prob_all <- function(time_intervals_list, rate_matrix_list, exact_time_rank = NULL) {

   mapply(
    FUN = function(time.intervals, rate.matrix) {
      probs.list <- lapply(time.intervals, FUN = "MatrixExp", mat = rate.matrix)
      
      if (!is.null(exact_time_rank)) {
        if (exact_time_rank < 1 || exact_time_rank > length(probs.list)) {
          stop("exact_time_rank must be between 1 and length(time.intervals)")
        }
        
        rate_nodiag <- rate.matrix
        diag(rate_nodiag) <- 0
        probs.list[[exact_time_rank]] <- probs.list[[exact_time_rank]] %*% rate_nodiag
      }
      
      probs.list
    },
    time.intervals = time_intervals_list,
    rate.matrix = rate_matrix_list,
    SIMPLIFY = FALSE
  )
}

#' Compute forward and backward probabilities for a hidden Markov model
#'
#' Computes forward probabilities, backward probabilities, and the observed-data
#' log-likelihood for a hidden Markov model.
#'
#' @param x Observed data sequence.
#' @param Pi List of transition probability matrices for intervals
#'   \(t_1 - t_2\), \(t_2 - t_3\), ..., \(t_{n-1} - t_n\).
#' @param delta Numeric vector giving the initial distribution of hidden states.
#' @param emission_matrix Matrix of emission probabilities. The \(i\)th row
#'   corresponds to the hidden state \(X(t)=i\), and the \(k\)th column gives
#'   \(P(O(t)=k \mid X(t)=i)\). Rows should sum to 1, and columns correspond to
#'   the possible observed states.
#'
#' @return A list with components:
#' \describe{
#'   \item{logalpha}{Log forward probabilities.}
#'   \item{logbeta}{Log backward probabilities.}
#'   \item{LL}{Observed-data log-likelihood.}
#' }
#'
forwardback <- function(x, Pi, delta, emission_matrix) {

  if (length(x) == 1) {
    LL <- as.numeric(log(sum(emission_matrix[, x] * delta)))
    return(list(LL = LL))
  }
  
  m <- nrow(Pi[[1]])
  n <- length(x)
  
  phi <- matrix(delta, nrow = 1)
  logalpha <- matrix(rep(NA_real_, m * n), nrow = n)
  lscale <- 0
  
  for (i in seq_len(n)) {
    if (i > 1) {
      phi <- phi %*% Pi[[i - 1]]
    }
    phi <- phi %*% diag(emission_matrix[, x[i]])
    sumphi <- as.numeric(sum(phi))
    phi <- phi / sumphi
    lscale <- lscale + log(sumphi)
    logalpha[i, ] <- as.numeric(log(phi) + lscale)
  }
  
  LL1 <- as.numeric(lscale)
  
  logbeta <- matrix(rep(NA_real_, m * n), nrow = n)
  logbeta[n, ] <- 0
  phi <- matrix(rep(1 / m, m), ncol = 1)
  lscale <- log(m)
  
  for (i in seq(n - 1, 1, -1)) {
    phi <- Pi[[i]] %*% diag(emission_matrix[, x[i + 1]]) %*% phi
    logbeta[i, ] <- as.numeric(log(phi) + lscale)
    sumphi <- as.numeric(sum(phi))
    phi <- phi / sumphi
    lscale <- lscale + log(sumphi)
  }
  
  list(
    logalpha = logalpha,
    logbeta = logbeta,
    LL = LL1
  )
}

#' Compute likelihood for one or more hidden Markov model observations
#'
#' Computes the likelihood across one or more observations under the hidden
#' Markov model.
#'
#' @param rates_list List of transition rate matrices.
#' @param init_list List of initial state distributions.
#' @param emission_list List of emission matrices.
#' @param obs_data_list List of observed state sequences.
#' @param obs_times_list List of observation times.
#' @param exact_time_rank Optional integer index of the interval with
#'   exact-event adjustment.
#'
#' @return The total likelihood across the given observations. 
likelihood <- function(rates_list, init_list, emission_list,
                       obs_data_list, obs_times_list,
                       exact_time_rank = NULL) {
 
  time_diffs_list <- lapply(obs_times_list, diff)
  
  transition_probabilities_list <- transition_prob_all(
    time_intervals_list = time_diffs_list,
    rate_matrix_list = rates_list,
    exact_time_rank = exact_time_rank
  )
  
  likelihood_forward_backward_list <- mapply(
    obs_data_list,
    transition_probabilities_list,
    init_list,
    emission_list,
    FUN = forwardback,
    SIMPLIFY = FALSE
  )
  
  LL <- vapply(likelihood_forward_backward_list, function(z) as.numeric(z$LL), numeric(1))
  sum(exp(LL))
}


#####################################################################
# 2. STANDARD TRIAL DESIGN INCIDENCE PROJECTION FUNCTIONS
#####################################################################
#' Compute control-arm incidence over an age interval under the standard design
#'
#' Computes control-arm incidence in a specified age interval under the
#' standard design.
#'
#' @param rate_matrix Transition rate matrix.
#' @param age1 Numeric value giving the start of the age interval.
#' @param age2 Numeric value giving the end of the age interval.
#' @param stage Character string indicating stage ("early" or "late").
#'
#' @return Interval incidence in the control arm.
get_control_incidence_standard <- function(rate_matrix, age1, age2, stage = c("early", "late")) {
  stage <- match.arg(stage)
  
  prob_mat_a1 <- MatrixExp(t = age1, mat = rate_matrix)
  prob_mat_a2 <- MatrixExp(t = age2, mat = rate_matrix)
  
  k <- nrow(rate_matrix)
  clinical_idx <- if (stage == "early") (k - 1) else k
  
  (prob_mat_a2[1, clinical_idx] - prob_mat_a1[1, clinical_idx]) /
    (1 - prob_mat_a1[1, k] - prob_mat_a1[1, k - 1])
}

#' Compute screen-arm incidence under the standard design
#'
#' Computes screen-arm incidence under the standard design for screen-detected,
#' interval-detected, or post-screen outcomes.
#'
#' @param sens_e Numeric value for early-stage sensitivity.
#' @param sens_l Numeric value for late-stage sensitivity.
#' @param screen_times Numeric vector of screening times.
#' @param age1 Numeric value giving the starting age.
#' @param rate_matrix Transition rate matrix.
#' @param stage Character string indicating stage ("early" or "late").
#' @param preclin Logical indicating whether the event is screen-detected.
#' @param post Logical indicating whether the event is post-screen clinical detection.
#' @param extended_followup Numeric value specifying additional follow-up time after the last screen.
#'
#' @return Cumulative incidence for the specified screen-arm detection type.    
get_screen_incidence_standard <- function(
   sens_e, sens_l,
    screen_times, age1, rate_matrix,
    stage = c("early", "late"),
    preclin = FALSE,
    post = FALSE,
    extended_followup = NULL
) {
  stage <- match.arg(stage)
  
  k <- nrow(rate_matrix)
  init_dist <- get_init(rate_matrix, age1)
  emission_dist <- get_emission(k, sens_e = sens_e, sens_l = sens_l)
  
  obs_data_list <- get_data(
    screen_times = screen_times,
    stage = stage,
    preclin = preclin,
    post = post,
    extended_followup = extended_followup
  )
  
  obs_times_list <- get_times(
    screen_times = screen_times,
    preclin = preclin,
    post = post,
    extended_followup = extended_followup
  )
  
  init_list     <- rep(list(init_dist), length(obs_data_list))
  emission_list <- rep(list(emission_dist), length(obs_data_list))
  rates_list    <- rep(list(rate_matrix), length(obs_data_list))
  
  likelihood(
    rates_list = rates_list,
    init_list = init_list,
    emission_list = emission_list,
    obs_data_list = obs_data_list,
    obs_times_list = obs_times_list
  )
}

#####################################################################
# 3. INTENDED EFFECT DESIGN INCIDENCE PROJECTION FUNCTIONS
#####################################################################
#' Compute true positivity rates across the screening period
#'
#' Computes early-stage, late-stage, and overall true positivity rates
#' (i.e., screen-positive probabilities) across the entire screening period.
#'
#' @param sens_e Numeric value for early-stage sensitivity.
#' @param sens_l Numeric value for late-stage sensitivity.
#' @param screen_times Numeric vector of screening times.
#' @param start_age Numeric value specifying the starting age.
#' @param rate_matrix Transition rate matrix.
#'
#' @return A list with components:
#' \describe{
#'   \item{positivity_early}{True positivity rate for early-stage disease.}
#'   \item{positivity_late}{True positivity rate for late-stage disease.}
#'   \item{positivity_overall}{Overall true positivity rate.}
#' }
get_positivity_rate <- function(sens_e, sens_l, screen_times, start_age, rate_matrix) {
    positivity_early <- get_screen_incidence_standard(
    sens_e = sens_e,
    sens_l = sens_l,
    screen_times = screen_times,
    age1 = start_age,
    rate_matrix = rate_matrix,
    stage = "early",
    preclin = TRUE
  )
  
  positivity_late <- get_screen_incidence_standard(
    sens_e = sens_e,
    sens_l = sens_l,
    screen_times = screen_times,
    age1 = start_age,
    rate_matrix = rate_matrix,
    stage = "late",
    preclin = TRUE
  )
  
  list(
    positivity_early = positivity_early,
    positivity_late = positivity_late,
    positivity_overall = positivity_early + positivity_late
  )
}

get_control_incidence_ie <- function(
  # --------------------------------------------------
    # Description:
    # Computes control-arm interval incidence under the intended-effect design,
    # conditional on ever-positive status.
    #
    # Inputs:
    # - sens_e: early-stage sensitivity
    # - sens_l: late-stage sensitivity
    # - screen_times: full vector of screen times
    # - start_age: starting age
    # - age1: start of interval
    # - age2: end of interval
    # - rate_matrix: transition rate matrix
    # - stage: "early" or "late"
    #
    # Outputs:
    # - Interval incidence in the control arm under the intended-effect design
    # --------------------------------------------------
    sens_e, sens_l,
    screen_times, start_age, age1, age2,
    rate_matrix,
    stage = c("early", "late")
) {
  stage <- match.arg(stage)
  clinical_code <- if (stage == "early") 4 else 5
  
  positivity_rate <- get_positivity_rate(
    sens_e = sens_e,
    sens_l = sens_l,
    screen_times = screen_times,
    start_age = start_age,
    rate_matrix = rate_matrix
  )$positivity_overall
  
  k <- nrow(rate_matrix)
  init_dist <- get_init(rate_matrix, start_age)
  emission_dist <- get_emission(k, sens_e = sens_e, sens_l = sens_l)
  
  t1 <- age1 - start_age
  t2 <- age2 - start_age
  
  if (t2 <= t1) {
    stop("age2 must be greater than age1.")
  }
  
  prob_first_positive_at_m_and_clinical_by_k <- function(m, k) {
    if (m < 1 || m > length(screen_times)) {
      stop("m is out of range.")
    }
    
    if (k< screen_times[m]) {
      return(0)
    }
    
    obs_times <- c(screen_times[1:m], k)
    
    first_positive_early <- c(rep(1, m - 1), 2, clinical_code)
    first_positive_late  <- c(rep(1, m - 1), 3, clinical_code)
    
    likelihood(
      rates_list = list(rate_matrix, rate_matrix),
      init_list = list(init_dist, init_dist),
      emission_list = list(emission_dist, emission_dist),
      obs_data_list = list(first_positive_early, first_positive_late),
      obs_times_list = list(obs_times, obs_times)
    )
  }
  
  valid_screens <- which(screen_times <= t1)
  
  if (length(valid_screens) == 0) {
    return(0)
  }
  
  interval_probability <- sum(
    sapply(valid_screens, function(m) {
      prob_first_positive_at_m_and_clinical_by_k(m, t2) -
        prob_first_positive_at_m_and_clinical_by_k(m, t1)
    })
  ) / positivity_rate
  
  interval_probability
}

get_screen_incidence_ie <- function(
    # --------------------------------------------------
    # Description:
    # Computes screen-detected incidence under the intended-effect design,
    # conditional on ever-positive status.
    #
    # Inputs:
    # - sens_e: early-stage sensitivity
    # - sens_l: late-stage sensitivity
    # - screen_times: vector of screen times
    # - start_age: starting age
    # - rate_matrix: transition rate matrix
    # - stage: "early" or "late"
    # - preclin: whether the event is screen-detected
    # - post: whether the event is postscreen clinical detection
    # - extended_followup: included for interface consistency
    # - positivity_rate: optional fixed overall positivity rate
    #
    # Outputs:
    # - Cumulative screen-detected incidence under the intended-effect design
    # --------------------------------------------------
    sens_e, sens_l,
    screen_times, start_age,
    rate_matrix,
    stage = c("early", "late"),
    preclin = FALSE,
    post = FALSE,
    extended_followup = NULL,
    positivity_rate = NULL
) {
  stage <- match.arg(stage)
  
  if (post) {
    return(0)
  }
  
  if (!preclin) {
    return(0)
  }
  
  if (is.null(positivity_rate)) {
    positivity_rate <- get_positivity_rate(
      sens_e = sens_e,
      sens_l = sens_l,
      screen_times = screen_times,
      start_age = start_age,
      rate_matrix = rate_matrix
    )$positivity_overall
  }
  
  get_screen_incidence_standard(
    sens_e = sens_e,
    sens_l = sens_l,
    screen_times = screen_times,
    age1 = start_age,
    rate_matrix = rate_matrix,
    stage = stage,
    preclin = TRUE
  ) / positivity_rate
}

#####################################################################
# 4. STAGE-SPEFIFIC TIME-DEPENDENT INCIDENCE PROJECTION FUNCTIONS BY TRIAL DESIGN AND ARM
#####################################################################

get_screen_arm_trial_results <- function(stage, design, params) {
  # --------------------------------------------------
  # Description:
  # Computes stage-specific screen-arm cancer incidence over the trial, separating
  # screen-detected, interval-detected, and screen+interval detected cancers.
  #
  # Inputs:
  # - stage: "early" or "late"
  # - design: "standard" or "IE"
  # - params: list with elements:
  #     start_age, numscreens, screen_int, num_followup_intervals,
  #     screen_times, rate_matrix, sens_e, sens_l, n_screen
  #
  # Outputs:
  # - list with:
  #     incidence:
  #       yearly, cumulative
  #     count:
  #       yearly, cumulative
  # - each contains:
  #     screen_detected
  #     interval_detected
  #     screen_arm_detected
  # --------------------------------------------------
  
  stage <- match.arg(stage, c("early", "late"))
  design <- match.arg(design, c("standard", "IE"))
  
  n_postscreen <- params$num_followup_intervals - 1
  n_intervals  <- params$numscreens + n_postscreen
  
  positivity_rate <- if (design == "IE") {
    get_positivity_rate(
      sens_e = params$sens_e,
      sens_l = params$sens_l,
      screen_times = params$screen_times,
      start_age = params$start_age,
      rate_matrix = params$rate_matrix
    )$positivity_overall
  } else {
    1
  }
  
  n_screen_effective <- params$n_screen * positivity_rate
  
  get_inc <- function(times, preclin = TRUE, post = FALSE, followup = NULL) {
    switch(
      design,
      standard = get_screen_incidence_standard(
        sens_e = params$sens_e,
        sens_l = params$sens_l,
        screen_times = times,
        age1 = params$start_age,
        rate_matrix = params$rate_matrix,
        stage = stage,
        preclin = preclin,
        post = post,
        extended_followup = followup
      ),
      IE = get_screen_incidence_ie(
        sens_e = params$sens_e,
        sens_l = params$sens_l,
        screen_times = times,
        start_age = params$start_age,
        rate_matrix = params$rate_matrix,
        stage = stage,
        preclin = preclin,
        post = post,
        extended_followup = followup,
        positivity_rate = positivity_rate
      )
    )
  }
  
  screen_detected_screen <- sapply(seq_len(params$numscreens), function(i) {
    get_inc(params$screen_times[1:i], preclin = TRUE)
  })
  
  cum_screen_detected <- c(
    screen_detected_screen,
    if (n_postscreen > 0) rep(tail(screen_detected_screen, 1), n_postscreen) else numeric(0)
  )
  
  cum_interval_detected <- switch(
    design,
    standard = {
      interval_screen <- sapply(seq_len(params$numscreens), function(i) {
        get_inc(params$screen_times[1:i], preclin = FALSE, post = FALSE, followup = params$screen_int)
      })
      
      interval_post <- if (n_postscreen > 0) {
        base_before_extra_followup <- if (params$numscreens > 1) {
          interval_screen[params$numscreens - 1]
        } else {
          0
        }
        
        sapply(seq_len(n_postscreen), function(j) {
          base_before_extra_followup +
            get_inc(
              params$screen_times,
              preclin = FALSE,
              post = TRUE,
              followup = (j + 1) * params$screen_int
            )
        })
      } else numeric(0)
      
      c(interval_screen, interval_post)
    },
    IE = rep(0, n_intervals)
  )
  
  cum_screen_arm_detected <- cum_screen_detected + cum_interval_detected
  
  yearly_screen_detected     <- diff(c(0, cum_screen_detected))
  yearly_interval_detected   <- diff(c(0, cum_interval_detected))
  yearly_screen_arm_detected <- yearly_screen_detected + yearly_interval_detected
  
  list(
    incidence = list(
      yearly = list(
        screen_detected     = yearly_screen_detected,
        interval_detected   = yearly_interval_detected,
        screen_arm_detected = yearly_screen_arm_detected
      ),
      cumulative = list(
        screen_detected     = cum_screen_detected,
        interval_detected   = cum_interval_detected,
        screen_arm_detected = cum_screen_arm_detected
      )
    ),
    count = list(
      yearly = list(
        screen_detected     = yearly_screen_detected * n_screen_effective,
        interval_detected   = yearly_interval_detected * n_screen_effective,
        screen_arm_detected = yearly_screen_arm_detected * n_screen_effective
      ),
      cumulative = list(
        screen_detected     = cum_screen_detected * n_screen_effective,
        interval_detected   = cum_interval_detected * n_screen_effective,
        screen_arm_detected = cum_screen_arm_detected * n_screen_effective
      )
    )
  )
}

get_control_arm_trial_results <- function(stage, design, params) {
  # --------------------------------------------------
  # Description:
  # Computes stage-specific control-arm cancer incidence over the trial under the
  # standard or intended-effect design.
  #
  # Inputs:
  # - stage: "early" or "late"
  # - design: "standard" or "IE"
  # - params: list with elements:
  #     start_age, numscreens, screen_int, num_followup_intervals,
  #     screen_times, rate_matrix, sens_e, sens_l, n_control
  #
  # Outputs:
  # - list with:
  #     incidence:
  #       yearly, cumulative
  #     count:
  #       yearly, cumulative
  # --------------------------------------------------
  stage <- match.arg(stage, c("early", "late"))
  design <- match.arg(design, c("standard", "IE"))
  
  n_intervals <- params$numscreens + params$num_followup_intervals - 1
  
  positivity_rate <- if (design == "IE") {
    get_positivity_rate(
      sens_e = params$sens_e,
      sens_l = params$sens_l,
      screen_times = params$screen_times,
      start_age = params$start_age,
      rate_matrix = params$rate_matrix
    )$positivity_overall
  } else {
    1
  }
  
  n_control_effective <- params$n_control * positivity_rate
  
  yearly_control_arm_detected <- mapply(
    function(age1, age2) {
      switch(
        design,
        standard = get_control_incidence_standard(
          rate_matrix = params$rate_matrix,
          age1 = age1,
          age2 = age2,
          stage = stage
        ),
        IE = get_control_incidence_ie(
          sens_e = params$sens_e,
          sens_l = params$sens_l,
          screen_times = params$screen_times,
          start_age = params$start_age,
          age1 = age1,
          age2 = age2,
          rate_matrix = params$rate_matrix,
          stage = stage
        )
      )
    },
    age1 = params$start_age + seq(0, by = params$screen_int, length.out = n_intervals),
    age2 = params$start_age + seq(params$screen_int, by = params$screen_int, length.out = n_intervals)
  )
  
  cum_control_arm_detected <- cumsum(yearly_control_arm_detected)
  
  list(
    incidence = list(
      yearly = list(
        control_arm_detected = yearly_control_arm_detected
      ),
      cumulative = list(
        control_arm_detected = cum_control_arm_detected
      )
    ),
    count = list(
      yearly = list(
        control_arm_detected = yearly_control_arm_detected * n_control_effective
      ),
      cumulative = list(
        control_arm_detected = cum_control_arm_detected * n_control_effective
      )
    )
  )
}


#####################################################################
# 5. MAIN FUNCTIONS
#####################################################################
#####################################################################
# 5.1. Main Function 1: Trial Projections (Both Standard and IE Designs)
# Generates results for both arms, all stages, and stage shift.
# Returns a list with all outcomes. 
#####################################################################

get_trial_results_by_design <- function(
    # --------------------------------------------------
    # Description:
    # Builds trial results for both standard and intended-effect designs, including
    # control-arm results, screen-arm results, and stage shift.
    #
    # Inputs:
    # - start_age: age at first screen
    # - numscreens: number of screens
    # - screen_int: interval between screens
    # - num_followup_intervals: number of follow-up intervals
    # - rate_matrix: natural history transition rate matrix
    # - sens_e: early-stage sensitivity
    # - sens_l: late-stage sensitivity
    # - n_control: control-arm size
    # - n_screen: screen-arm size
    #
    # Outputs:
    # - list with:
    #     params: input values
    #     timeline: time intervals with ages and period labels
    #     control: results by design ("standard", "IE") and stage ("early","late","overall")
    #     screen: results by design ("standard", "IE") and stage ("early","late","overall")
    #     stage_shift: yearly and cumulative late-stage stage shift
    # --------------------------------------------------
    start_age,
    numscreens,
    screen_int,
    num_followup_intervals,
    rate_matrix,
    sens_e,
    sens_l,
    n_control,
    n_screen
) {
  
  params <- list(
    start_age = start_age,
    numscreens = numscreens,
    screen_int = screen_int,
    num_followup_intervals = num_followup_intervals,
    screen_times = seq(0, (numscreens - 1) * screen_int, by = screen_int),
    rate_matrix = rate_matrix,
    sens_e = sens_e,
    sens_l = sens_l,
    n_control = n_control,
    n_screen = n_screen
  )
  
  n_intervals <- numscreens + num_followup_intervals - 1
  
  timeline <- data.frame(
    time_since_randomization = seq_len(n_intervals),
    age1 = start_age + (seq_len(n_intervals) - 1) * screen_int,
    age2 = start_age + seq_len(n_intervals) * screen_int,
    period = c(
      rep("screen", numscreens),
      rep("postscreen", num_followup_intervals - 1)
    )
  )
  
  control <- list(
    standard = list(
      early = get_control_arm_trial_results("early", "standard", params),
      late  = get_control_arm_trial_results("late",  "standard", params)
    ),
    IE = list(
      early = get_control_arm_trial_results("early", "IE", params),
      late  = get_control_arm_trial_results("late",  "IE", params)
    )
  )
  
  add_control_overall <- function(x) {
    x$overall <- list(
      incidence = list(
        yearly = list(
          control_arm_detected =
            x$early$incidence$yearly$control_arm_detected +
            x$late$incidence$yearly$control_arm_detected
        ),
        cumulative = list(
          control_arm_detected =
            x$early$incidence$cumulative$control_arm_detected +
            x$late$incidence$cumulative$control_arm_detected
        )
      ),
      count = list(
        yearly = list(
          control_arm_detected =
            x$early$count$yearly$control_arm_detected +
            x$late$count$yearly$control_arm_detected
        ),
        cumulative = list(
          control_arm_detected =
            x$early$count$cumulative$control_arm_detected +
            x$late$count$cumulative$control_arm_detected
        )
      )
    )
    x
  }
  
  control$standard <- add_control_overall(control$standard)
  control$IE <- add_control_overall(control$IE)
  
  screen <- list(
    standard = list(
      early = get_screen_arm_trial_results("early", "standard", params),
      late  = get_screen_arm_trial_results("late",  "standard", params)
    ),
    IE = list(
      early = get_screen_arm_trial_results("early", "IE", params),
      late  = get_screen_arm_trial_results("late",  "IE", params)
    )
  )
  
  add_screen_overall <- function(x) {
    x$overall <- list(
      incidence = list(
        yearly = list(
          screen_detected =
            x$early$incidence$yearly$screen_detected +
            x$late$incidence$yearly$screen_detected,
          interval_detected =
            x$early$incidence$yearly$interval_detected +
            x$late$incidence$yearly$interval_detected,
          screen_arm_detected =
            x$early$incidence$yearly$screen_arm_detected +
            x$late$incidence$yearly$screen_arm_detected
        ),
        cumulative = list(
          screen_detected =
            x$early$incidence$cumulative$screen_detected +
            x$late$incidence$cumulative$screen_detected,
          interval_detected =
            x$early$incidence$cumulative$interval_detected +
            x$late$incidence$cumulative$interval_detected,
          screen_arm_detected =
            x$early$incidence$cumulative$screen_arm_detected +
            x$late$incidence$cumulative$screen_arm_detected
        )
      ),
      count = list(
        yearly = list(
          screen_detected =
            x$early$count$yearly$screen_detected +
            x$late$count$yearly$screen_detected,
          interval_detected =
            x$early$count$yearly$interval_detected +
            x$late$count$yearly$interval_detected,
          screen_arm_detected =
            x$early$count$yearly$screen_arm_detected +
            x$late$count$yearly$screen_arm_detected
        ),
        cumulative = list(
          screen_detected =
            x$early$count$cumulative$screen_detected +
            x$late$count$cumulative$screen_detected,
          interval_detected =
            x$early$count$cumulative$interval_detected +
            x$late$count$cumulative$interval_detected,
          screen_arm_detected =
            x$early$count$cumulative$screen_arm_detected +
            x$late$count$cumulative$screen_arm_detected
        )
      )
    )
    x
  }
  
  screen$standard <- add_screen_overall(screen$standard)
  screen$IE <- add_screen_overall(screen$IE)
  
  get_stage_shift <- function(design) {
    control_late_yearly <- control[[design]]$late$incidence$yearly$control_arm_detected
    screen_late_yearly  <- screen[[design]]$late$incidence$yearly$screen_arm_detected
    
    control_late_cumulative <- control[[design]]$late$incidence$cumulative$control_arm_detected
    screen_late_cumulative  <- screen[[design]]$late$incidence$cumulative$screen_arm_detected
    
    list(
      yearly = ifelse(
        control_late_yearly > 0,
        100 * (control_late_yearly - screen_late_yearly) / control_late_yearly,
        NA_real_
      ),
      cumulative = ifelse(
        control_late_cumulative > 0,
        100 * (control_late_cumulative - screen_late_cumulative) / control_late_cumulative,
        NA_real_
      )
    )
  }
  
  stage_shift <- list(
    standard = get_stage_shift("standard"),
    IE = get_stage_shift("IE")
  )
  
  round_nested <- function(x, digits) {
    if (is.numeric(x)) {
      round(x, digits)
    } else if (is.list(x)) {
      lapply(x, round_nested, digits = digits)
    } else {
      x
    }
  }
  
  round_count_only <- function(x) {
    if (!is.list(x)) return(x)
    
    for (nm in names(x)) {
      if (nm == "count") {
        x[[nm]] <- round_nested(x[[nm]], 0)
      } else if (is.list(x[[nm]])) {
        x[[nm]] <- round_count_only(x[[nm]])
      }
    }
    x
  }
  
  ## round ONLY counts
  control <- round_count_only(control)
  screen  <- round_count_only(screen)
  
  ## round stage shift
  stage_shift <- round_nested(stage_shift, 0)
  list(
    params = params,
    timeline = timeline,
    control = control,
    screen = screen,
    stage_shift = stage_shift
  )
}

#####################################################################
# 5.2. Main Function: User-Friendly Output Tables
# Generates formatted tables from the results list produced by
# Main Function 1, with options for trial design, stage and outcome type selection.
#####################################################################

stage_shift_by_design <- function(
    # --------------------------------------------------
    # Description:
    # Creates a table of late-stage results and stage shift for one design, with
    # optional early-stage and overall results.
    #
    # Inputs:
    # - results: output from get_trial_results_by_design()
    # - design: "standard" or "IE"
    # - yearly_or_cumulative: "yearly" or "cumulative"
    # - incidence_or_count: "incidence" or "count"
    # - early_stage: TRUE/FALSE (include early-stage results)
    # - overall: TRUE/FALSE (include overall results)
    #
    # Outputs:
    # - data frame with:
    #     time_since_randomization, age1, age2, period
    #     late-stage control and screen results
    #     optional early-stage results
    #     optional overall results
    #     stage_shift
    # --------------------------------------------------
    results,
    design = c("standard", "IE"),
    yearly_or_cumulative = c("yearly", "cumulative"),
    incidence_or_count = c("incidence", "count"),
    early_stage = FALSE,
    overall = FALSE
) {
  
  design <- match.arg(design)
  yearly_or_cumulative <- match.arg(yearly_or_cumulative)
  incidence_or_count <- match.arg(incidence_or_count)
  
  get_stage_cols <- function(stage_name, suffix = stage_name) {
    data.frame(
      setNames(
        list(results$control[[design]][[stage_name]][[incidence_or_count]][[yearly_or_cumulative]]$control_arm_detected),
        paste0("control_arm_detected_", suffix)
      ),
      setNames(
        list(results$screen[[design]][[stage_name]][[incidence_or_count]][[yearly_or_cumulative]]$screen_detected),
        paste0("screen_detected_", suffix)
      ),
      setNames(
        list(results$screen[[design]][[stage_name]][[incidence_or_count]][[yearly_or_cumulative]]$interval_detected),
        paste0("interval_detected_", suffix)
      ),
      setNames(
        list(results$screen[[design]][[stage_name]][[incidence_or_count]][[yearly_or_cumulative]]$screen_arm_detected),
        paste0("screen_arm_detected_", suffix)
      )
    )
  }
  
  out <- data.frame(
    time_since_randomization = results$timeline$time_since_randomization,
    age1 = results$timeline$age1,
    age2 = results$timeline$age2,
    period = results$timeline$period
  )
  
  out <- cbind(out, get_stage_cols("late", "late"))
  
  if (isTRUE(early_stage)) {
    out <- cbind(out, get_stage_cols("early", "early"))
  }
  
  if (isTRUE(overall)) {
    out <- cbind(out, get_stage_cols("overall", "overall"))
  }
  
  out$stage_shift <- results$stage_shift[[design]][[yearly_or_cumulative]]
  
  out
}