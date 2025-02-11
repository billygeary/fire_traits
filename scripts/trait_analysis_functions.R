# Functions

plot_pcoa <- function(data.in) {
  # Calculate convex hulls for each group
  hulls <- data.in %>%
    group_by(Cluster) %>%  # Use tidy evaluation for dynamic grouping
    slice(chull(NMDS1, NMDS2))           # Select points forming the convex hull
  
  # Create the PCA plot
  pcoa_plot <- ggplot(data.in) + 
    geom_point(aes(x = NMDS1, y = NMDS2, colour = Cluster), size = 4) + 
    geom_text(aes(x = NMDS1, y = NMDS2, label = Species), nudge_y = -0.01, size = 3) +
    geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, 
                                   group = Cluster, 
                                   fill = Cluster), alpha = 0.3) +
    theme_bw() +
    xlab("NMDS1") + 
    ylab("NMDS2") #+ facet_wrap(~Cluster, scales="free")
  
  return(pcoa_plot)
}


plot_pcoa_with_vectors <- function(data.in, vectors) {
  # Calculate convex hulls for each group
  hulls <- data.in %>%
    group_by(Cluster) %>% 
    slice(chull(NMDS1, NMDS2))
  
  # Create NMDS plot with vectors
  pcoa_plot <- ggplot(data.in) + 
    geom_point(aes(x = NMDS1, y = NMDS2, colour = Cluster), size = 4) + 
    geom_segment(data = vectors, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                 arrow = arrow(length = unit(0.2, "cm")), color = "black") + 
    geom_text(data = vectors, aes(x = NMDS1, y = NMDS2, label = rownames(vectors)),
              color = "black", hjust = 0, vjust = -0.5, position=position_dodge2(width=1)) + 
    coord_cartesian(ylim=c(-0.6,0.6), xlim=c(-0.5, 0.9)) +
    theme_bw() +
    xlab("NMDS1") + 
    ylab("NMDS2")
  
  return(pcoa_plot)
}

plot_pcoa_aquatic <- function(data.in) {
  # Calculate convex hulls for each group
  hulls <- data.in %>%
    group_by(Cluster) %>%  # Use tidy evaluation for dynamic grouping
    slice(chull(NMDS1, NMDS2))           # Select points forming the convex hull
  
  # Create the PCA plot
  pcoa_plot <- ggplot(data.in) + 
    geom_point(aes(x = NMDS1, y = NMDS2, colour = Cluster), size = 4) + 
    geom_text(aes(x = NMDS1, y = NMDS2, label = Species), nudge_y = -0.01) +
    geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, 
                                   group = Cluster, 
                                   fill = Cluster), alpha = 0.3) +
    theme_bw() + #facet_wrap(~Cluster, scales="free", nrow=2) +
    xlab("NMDS1") + 
    ylab("NMDS2") + facet_wrap(~aquatic, scales = "free")
  
  return(pcoa_plot)
}

# Extract partial dependence data for each predictor
get_pdp_data <- function(var, model) {
  pd <- plot.gbm(model, i.var = var, return.grid = TRUE)
  pd$variable <- var  # Add variable name
  colnames(pd) <- c("x", "y", "variable")
  # Convert categorical variables to factors for proper ggplot treatment
  if (var %in% factor_vars) {
    pd$x <- as.factor(pd$x)
  }
  return(pd)
}

# Create ggplot for each variable
create_pdp_plots = function(df.in, factor_vars){
  df<-as.data.frame(df.in)
  if (df$variable[1] %in% factor_vars) {
    # Bar plot for categorical variables
    plot <- ggplot(df, aes(x = x, y = y, fill = x)) +
      geom_col(show.legend = FALSE) +
      labs(title = df$variable[1], x = df$variable[1], y = "Partial Dependence") +
      theme_minimal()
  } else {
    # Line plot for continuous variables
    plot <- ggplot(df, aes(x = x, y = y)) +
      geom_line(color = "blue") +
      labs(title = df$variable[1], x = df$variable[1], y = "Partial Dependence") +
      theme_minimal()
  }
  return(plot)
}

# Function to get predicted effects for a given predictor
get_pred_data <- function(var, model) {
  pred <- ggeffects::ggpredict(model, terms = var)  # Get predicted effects
  pred$variable <- var  # Store variable name
  return(pred)
}


create_pred_plots = function(df.in, factor_vars){
  df<-as.data.frame(df.in)
  if (df$variable[1] %in% factor_vars) {
    # Bar plot for categorical variables
   plot <- ggplot(df, aes(x = x, y = predicted, ymin=conf.low,ymax=conf.high, fill = x)) +
      geom_pointrange(show.legend = FALSE) +
      labs(title = df$variable[1], x = df$variable[1], y = "Predicted Response") +
      theme_minimal()
  } else {
    # Line plot for continuous variables
  plot <- ggplot(df, aes(x = x, y = predicted)) +
      geom_line(color = "black") +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +  # Add confidence interval
      labs(title = df$variable[1], x = df$variable[1], y = "Predicted Response") +
      theme_minimal()
  }
  return(plot)
}



