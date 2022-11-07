rm(list = ls()) # emtpy workspace

source("BSO_main_function.R")

# load Holzinger & Swineford data set
data("HS.ability.data", package = "OpenMx")
HS_data <- HS.ability.data
item_names <- colnames(HS_data[,7:30])
HS_data <- as.data.frame(apply(HS_data[, item_names], 2, scale))


# This is a small working example which will run several minutes (< 10 min.)
# For more realistic settings of the hyperparameters please see text
max_iter <- 10     		 # number of iterations without improvement after which to abort search
n_bees <- 200     		 # number of n_bees
top_best <- 40      	 # how many solutions to use for onlookers
percent_scouts <- .25  # percentage scouts
min_nest_fac <- 1    	 # how many nested factors minimum
max_nest_fac <- 5    	 # how many nested factors maximum
depletion <- 5 			   # update best solutions after how many iterations
cluster_mode <- FALSE  # Are you running on cluster? If so, do not repeat plotting every iteration


for (seed in 1) {

  res.name <- paste0("BSO_s" , seed , "_b" , n_bees, "_s", percent_scouts*n_bees, "_t", top_best, Sys.Date())
  summaryfile <- paste0(res.name, ".csv")
  summaryfile_fin = paste0(res.name, "_final.csv")
  
  conv_plot <- 
    BSO(item_names = item_names, 
        data = HS_data,
        max_iter = max_iter, 
        n_bees = n_bees,
        n_start_bees = n_bees,
        percent_scouts = percent_scouts, 
        top_best = top_best, 
        min_nest_fac = min_nest_fac, 
        max_nest_fac = max_nest_fac, 
        depletion = depletion, 
        summaryfile = summaryfile,
        summaryfile_fin = summaryfile_fin,
        seed = seed,
        verbose = FALSE,
        parallel = TRUE,
        nCores = 10,
        plot_nectar = TRUE,
        fit_crit = c("cfi", "rmsea", "min_loading", "n_items"),
        logistic_weights = list(c(d = 0.95, a = 50),
                                c(d = 0.06, a = -50),
                                c(d = 0.3, a = 15),
                                c(d = 0.95, a = 30)),
        nu_weights = c(2,2,1,2),
        nu_min = 0.0001,
        balance_n_fac = FALSE,
        n_start_items = 24,
        bounds = "pos.var",
        cluster_mode = cluster_mode,
        plot_list = list(xlim = c(0, max_iter*10), 
                         ylim = c(0, sum(c(2,2,1,2))),
                         ylab = "Overall Nectar Value",
                         xlab = "Iteration",
                         jitter_width = 0.5,
                         alpha = 0.2,
                         size = 1)
    )
  
  if (!is.null(conv_plot)){
    suppressWarnings(expr = {
      factor_cols <- hue_pal()(5)
      if(!cluster_mode) plot(conv_plot)
      ggsave(plot = conv_plot, filename = paste0(res.name,".pdf"), device = "pdf")
    })
  }
}
