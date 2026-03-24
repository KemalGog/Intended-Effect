#Power functions

get_p_trad <- function(p_c, p_s, R = 1) {
  (p_c + R * p_s) / (1 + R)
}

get_p_ie <- function(p_c_ie, p_s_ie, R = 1) {
  (p_c_ie + R * p_s_ie) / (1 + R)
}

risk_diff_trad <- function(p_c, p_s) {
  p_c - p_s
}

risk_diff_ie <- function(p_c_ie, p_s_ie) {
  p_c_ie - p_s_ie
}

efficiency_ratio_ie_vs_trad_from_probs <- function(pi_pos, p_c, p_s, p_c_ie, p_s_ie, R = 1) {
  p <- get_p_trad(p_c, p_s, R)
  p_ie <- get_p_ie(p_c_ie, p_s_ie, R)
  (1 / pi_pos) * (p_ie * (1 - p_ie)) / (p * (1 - p))
}

power_trad <- function(n_total, p_c, p_s, alpha = 0.025, R = 1) {
  n_c <- n_total / (1 + R)
  n_s <- n_total - n_c
  p <- get_p_trad(p_c = p_c, p_s = p_s, R = R)
  z_alpha <- qnorm(1 - alpha)
  se <- sqrt(p * (1 - p) * (1 / n_c + 1 / n_s))
  z_effect <- (p_c - p_s) / se
  1 - pnorm(z_alpha - z_effect)
}

n_trad_required <- function(power = 0.80, p_c, p_s, alpha = 0.025, R = 1) {
  p <- get_p_trad(p_c = p_c, p_s = p_s, R = R)
  z_alpha <- qnorm(1 - alpha)
  z_beta <- qnorm(power)
  ((1 + R)^2 * (z_alpha + z_beta)^2 * p * (1 - p)) /
    (R * (p_c - p_s)^2)
}

power_ie <- function(n_total, pi_pos, p_c_ie, p_s_ie, alpha = 0.025, R = 1) {
  n_c <- n_total / (1 + R)
  n_s <- n_total - n_c
  p_ie <- get_p_ie(p_c_ie = p_c_ie, p_s_ie = p_s_ie, R = R)
  z_alpha <- qnorm(1 - alpha)
  se <- sqrt((1 / pi_pos) * p_ie * (1 - p_ie) * (1 / n_c + 1 / n_s))
  z_effect <- (p_c_ie - p_s_ie) / se
  1 - pnorm(z_alpha - z_effect)
}

n_ie_required <- function(power = 0.80, pi_pos, p_c_ie, p_s_ie, alpha = 0.025, R = 1) {
  p_ie <- get_p_ie(p_c_ie = p_c_ie, p_s_ie = p_s_ie, R = R)
  z_alpha <- qnorm(1 - alpha)
  z_beta <- qnorm(power)
  (1 / pi_pos) *
    ((1 + R)^2 * (z_alpha + z_beta)^2 * p_ie * (1 - p_ie)) /
    (R * (p_c_ie - p_s_ie)^2)
}

#FIGURE HELPER FUNCTION
#--------------------------------------------------
# Helper: choose nice common y-axis for incidence/count panels
#--------------------------------------------------
get_detection_ymax <- function(df, show_early = FALSE, show_overall = FALSE) {
  
  vars <- c("control_arm_detected_late", "screen_arm_detected_late")
  
  if (show_early) {
    vars <- c(vars, "control_arm_detected_early", "screen_arm_detected_early")
  }
  
  if (show_overall) {
    vars <- c(vars, "control_arm_detected_overall", "screen_arm_detected_overall")
  }
  
  ymax <- max(unlist(df[vars]), na.rm = TRUE)
  return(ymax)
}

get_achieved_power <- function(n_required_fun, N_current, interval = c(0.001, 0.999)) {
  
  f <- function(p) n_required_fun(p) - N_current
  
  f_low  <- f(interval[1])
  f_high <- f(interval[2])
  
  if (f_low < 0 && f_high < 0) {
    return(">0.999")
  } else if (f_low > 0 && f_high > 0) {
    return("<0.001")
  } else {
    return(round(uniroot(f, interval = interval)$root, 3))
  }
}

#--------------------------------------------------
# Helper: common theme
#--------------------------------------------------
library(ggplot2)
library(patchwork)
library(scales)

get_detection_ymax <- function(df, show_early = FALSE, show_overall = FALSE) {
  vars <- c("control_arm_detected_late", "screen_arm_detected_late")
  
  if (show_early) {
    vars <- c(vars, "control_arm_detected_early", "screen_arm_detected_early")
  }
  
  if (show_overall) {
    vars <- c(vars, "control_arm_detected_overall", "screen_arm_detected_overall")
  }
  
  max(unlist(df[vars]), na.rm = TRUE)
}

theme_trial_panels <- function() {
  theme_minimal(base_size = 13) +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      legend.position = "top",
      legend.title = element_blank(),
      
      # remove grid
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # add axis lines + ticks
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.2, "cm")
    )
}

make_detection_plot <- function(df, y_control, y_screen, y_lab,
                                ymax_det, count_or_incidence) {
  
  ggplot(df, aes(x = time_since_randomization)) +
    geom_line(aes(y = .data[[y_control]], color = "Control"), linewidth = 1.2) +
    geom_line(aes(y = .data[[y_screen]], color = "Screen"), linewidth = 1.2) +
    labs(
      x = "Time since randomization (years)",
      y = y_lab,
      color = NULL
    ) +
    scale_x_continuous(
      breaks = seq(
        floor(min(df$time_since_randomization, na.rm = TRUE)),
        ceiling(max(df$time_since_randomization, na.rm = TRUE)),
        by = 1
      )
    ) +
    coord_cartesian(ylim = c(0, ymax_det)) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.02)),
      labels = if (count_or_incidence == "incidence") {
        scales::label_number(accuracy = 0.001)
      } else {
        waiver()
      }
    ) +
    theme_trial_panels()
}

make_stage_shift_plot <- function(df, cumulative_or_interval) {
  
  ymin_shift <- floor(min(df$stage_shift, na.rm = TRUE) / 10) * 10
  
  ggplot(df, aes(x = time_since_randomization, y = stage_shift)) +
    geom_line(linewidth = 1.2) +
    labs(
      x = "Time since randomization (years)",
      y = "Stage shift (%)"
    ) +
    scale_x_continuous(
      breaks = seq(
        floor(min(df$time_since_randomization, na.rm = TRUE)),
        ceiling(max(df$time_since_randomization, na.rm = TRUE)),
        by = 1
      )
    ) +
    scale_y_continuous(
      limits = c(ymin_shift, 100),
      breaks = seq(ymin_shift, 100, by = 10)
    ) +
    theme_trial_panels() +
    theme(legend.position = "none")
}

plot_trial_panels <- function(trial_table_all,
                              show_early = FALSE,
                              show_overall = FALSE,
                              cumulative_or_interval = "cumulative",
                              count_or_incidence = "count",
                              design_name = "Standard Trial Design") {
  
  ymax_det <- get_detection_ymax(
    trial_table_all,
    show_early = show_early,
    show_overall = show_overall
  )
  
  y_lab_late <- ifelse(
    count_or_incidence == "count",
    "Number of late-stage cancers detected",
    "Late-stage incidence"
  )
  
  y_lab_early <- ifelse(
    count_or_incidence == "count",
    "Number of early-stage cancers detected",
    "Early-stage incidence"
  )
  
  y_lab_overall <- ifelse(
    count_or_incidence == "count",
    "Number of cancers detected",
    "Overall incidence"
  )
  
  p_late <- make_detection_plot(
    df = trial_table_all,
    y_control = "control_arm_detected_late",
    y_screen  = "screen_arm_detected_late",
    y_lab = y_lab_late,
    ymax_det = ymax_det,
    count_or_incidence = count_or_incidence
  )
  
  p_shift <- make_stage_shift_plot(
    df = trial_table_all,
    cumulative_or_interval = cumulative_or_interval
  )
  
  plot_list <- list(p_late, p_shift)
  
  if (show_early) {
    p_early <- make_detection_plot(
      df = trial_table_all,
      y_control = "control_arm_detected_early",
      y_screen  = "screen_arm_detected_early",
      y_lab = y_lab_early,
      ymax_det = ymax_det,
      count_or_incidence = count_or_incidence
    ) + theme(legend.position = "none")
    
    plot_list <- c(plot_list, list(p_early))
  }
  
  if (show_overall) {
    p_overall <- make_detection_plot(
      df = trial_table_all,
      y_control = "control_arm_detected_overall",
      y_screen  = "screen_arm_detected_overall",
      y_lab = y_lab_overall,
      ymax_det = ymax_det,
      count_or_incidence = count_or_incidence
    ) + theme(legend.position = "none")
    
    plot_list <- c(plot_list, list(p_overall))
  }
  
  overall_title <- paste0(
    time_label(cumulative_or_interval),
    " Outcomes Under ",
    design_name
  )
  
  wrap_plots(plotlist = plot_list, ncol = 2) +
    plot_annotation(
      title = overall_title,
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
      )
    )
}