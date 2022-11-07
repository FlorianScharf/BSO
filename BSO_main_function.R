# Load all required packages
library(psych)
library(lavaan)
library(parallel)
library(GPArotation)
library(OpenMx)
library(semTools)
library(ggplot2)
library(scales)
# Load necessary helper functions
source("tool_functions.R")
source("fit_function.R")



########################## Iterations ################################################
BSO = function(item_names, # a vector of item names (length must equal the length of item assignments)
               data, # a data frame containing the data set
               max_iter, # how many iterations should BSO try to find a better solution than the current optimum?
               n_bees, # how many bees to use
               n_start_bees = n_bees, # how many random starts? as many as bees if no value provided
               n_start_items = length(item_names), # how many items should be in the start model
               percent_scouts, # proportion of scouts among all bees
               top_best, # how many of the best solutions 
               min_nest_fac, # minimum number of nested factors 
               max_nest_fac, # maximum number of nested factors
               depletion, # how long should BSO try to improve one specific solution before giving up on it
               summaryfile, # name of the detailed summary file
               summaryfile_fin, # name of a shortened summary file with best model only
               scouts = round(prop.table(top_best:1) * n_bees * percent_scouts), #scouts per top solution
               onlookers = round(prop.table(top_best:1) * n_bees * (1 - percent_scouts)),  #onlookers per top solution
               fit_crit = c("cfi", "rmsea", "min_omega3", "min_loading"), # a vector of fit criteria to be extracted
               # list (length equal to fit_crit) with weights for the logistic transformation (see section The Optimization Function)
               logistic_weights = list(c(d = 0.9, a = 70),
                                       c(d = 0.06, a = -70),
                                       c(d = 0.4, a = 70),
                                       c(d = 0.33, a = 70)),
               nu_weights = c(1,1,1,1), # a vector of weights for the aggregation function (see section The Optimization Function)
               nu_min = 0, # threshold value, if the nu value of any fit criterion is below this threshold, nu
               # will be returned as zero to prevent BSO from optimizing a subset of the fit criteria while neglecting
               # one criterion completely
               verbose = FALSE, # return information in console (incompatible with parallel mode!)
               debug_fit_mode = FALSE,  # return fit object for debugging
               parallel = TRUE, # spread computations across CPU cores?
               nCores = parallel::detectCores() - 2, # Number of cores for parallel mode
               seed = 1, # random seed to be used 
               fun = TRUE, # see for yourself...
               plot_nectar = TRUE, # create a plot of the nectar value after during computation
               # some aesthetic parameters for the nectar plot
               plot_list = list(xlim = c(0,max_iter*2), 
                                ylim = c(0,sum(nu_weights)),
                                ylab = "Overall Nectar Value",
                                xlab = "Iteration",
                                jitter_width = 0.5,
                                alpha = 0.2,
                                size = 1),
               cluster_mode = FALSE, # for use on computing cluster: deactivates the plot to not spam a PDF-plot-file...
               balance_n_fac = FALSE, # # if TRUE, the start models are sampled such that the number of factors is uniformly distributed
               ...){ # further arguments to be passed to lavaan 
  
  # see for yourself
  if(fun){
    cat(readLines("fun1.txt", warn = F), sep = "\n")
  }
  
 # Initialize nectar plot
    if (plot_nectar){

    conv_plot <- ggplot() + 
      xlim(plot_list$xlim[1], plot_list$xlim[2]) +
      ylim(plot_list$ylim[1], plot_list$ylim[2]) +
      xlab(plot_list$xlab) +
      ylab(plot_list$ylab) 
    if(!cluster_mode) plot(conv_plot)
      }
  
  
  # check if the number of bees is consistent with the computed number of
  # scouts and onlookers
  if (sum(scouts + onlookers) != n_bees) {
    n_bees <- sum(scouts + onlookers)
    warning(paste0("The current combination of n_bees, percent_scouts and top_best does not work out, I changed n_bees to "), n_bees)
  }
  
  # in parallel mode initialize PSOCKET-Cores (for windows compatibility)
  if (parallel) {
    myCL <- makePSOCKcluster(names = nCores, outfile = ifelse(verbose, paste0("parallel_verbose",Sys.time(),".txt"), ""))
    clusterSetRNGStream(cl = myCL, iseed = seed)
    clusterExport(cl = myCL, varlist = ls(envir = environment()), envir = environment())
    clusterExport(cl = myCL, varlist = ls(envir = .GlobalEnv), envir = .GlobalEnv)
    clusterEvalQ(cl = myCL, expr = {
      library(psych)
      library(lavaan)
      library(semTools)
      library(ggplot2)
      library(scales)
    })
   
  } else {set.seed(seed)
      }
  
  n_items <- length(item_names) #number of items in total
  
  ########################### Initial Scouts ##############################
  # Initial scouts will randomly select models
  
 if (!parallel){
   solutions <- sapply(1:n_start_bees, function(i_scout){ 
     if (verbose) print(i_scout)
     initialize.bees(item_names = item_names,
                     data = data,
                     max_iter = max_iter,
                     n_bees = n_bees,
                     n_start_items = n_start_items,
                     fit_crit = fit_crit,
                     logistic_weights = logistic_weights,
                     nu_weights = nu_weights,
                     nu_min = nu_min,
                     balance_n_fac = balance_n_fac,
                     debug_fit_mode = debug_fit_mode,
                     verbose = verbose,
                     ...)
   }, simplify = "matrix")
 } else if (parallel){
  solutions <- parSapply(myCL, 1:n_start_bees, function(i_scout){ 
    if (verbose) print(i_scout)
    initialize.bees(item_names = item_names,
                    data = data,
                    max_iter = max_iter,
                    n_bees = n_bees,
                    n_start_items = n_start_items,
                    fit_crit = fit_crit,
                    logistic_weights = logistic_weights,
                    nu_weights = nu_weights,
                    nu_min = nu_min,
                    balance_n_fac = balance_n_fac,
                    debug_fit_mode = debug_fit_mode,
                    verbose = verbose,
                    ...
                    )
   }, simplify = "matrix")
 }

  # Transpose solutions object for convenience purposes
  solutions <- data.frame(t(solutions))
  # Sort solutions by nu so that the top_best are on top
  solutions <- solutions[order(solutions$fit_overall, decreasing = TRUE),]
  
  ################################################ Scouts & Onlookers ##############################################

  #onlookers will modify the best overall models
  iter <- 1    #overall iteration counter
  counter <- 1    #iteration counter that will be reset when a improved solution is found
  best_fit <- solutions$fit_overall[1] #used to store best fit of all iterations
  best_solution <- solutions[1,]


  while (counter <= max_iter){ #while resetable count <= maximum iteration number
    start <- Sys.time() # measure run time
    
    tmp_solutions <- solutions[1:top_best,] # solutions to be investigated further
    fresh <- solutions$age <= depletion # check for depletion
    solutions$age[fresh][1:top_best] <- solutions$age[fresh][1:top_best] + 1 # increase "age" of best solutions
  
    # Prepare a data frame with the solutions to be modified in this iteration
    allBees <- data.frame(i_bee_ind = 1:n_bees)
    allBees$solution_focus = c(rep(c(1:top_best), times = scouts), rep(c(1:top_best), times = onlookers)) 
    allBees$role = c(rep("scout", sum(scouts)), rep("onlooker", sum(onlookers)))
    allBees <- data.frame(allBees, tmp_solutions[allBees$solution_focus,item_names]) 
   
    # go through the models in parallel or serial mode
    if (!parallel){
      new_solutions <- sapply(1:nrow(allBees), function(i_bee){
        
        if (verbose) print(allBees[i_bee,])
        
        item_assignment_orig <- allBees[i_bee, item_names]
      
        # change model depending on the role of the current bee
        if (allBees$role[i_bee] == "scout"){
          item_assignment <- scout.bee(item_assignment = item_assignment_orig,
                                       max_nest_fac = max_nest_fac)
          
          
        } else if (allBees$role[i_bee] == "onlooker") {
          
          item_assignment <- onlooker.bee(item_assignment = item_assignment_orig, 
                                          max_nest_fac = max_nest_fac)
          
        } else {stop("Fatal error: Bee swarm out of control.")} # if an unknown role occurs...
        
        #checks if this model has been tested before
        compare_item_assignment <- apply(solutions[, item_names], 1, function(i_solution) identical(item_assignment, i_solution))

        if (any(compare_item_assignment)){ # if you saw this item configuration before
          if (verbose) print ("Empty Flower")
          new_sol <- solutions[max(which(compare_item_assignment)),]
          # the new age is the maximal age of this solution + 1
          # (because solutions is sorted by fit)
          new_sol$age <- max(solutions$age[which(compare_item_assignment)]) + 1 
        } else { # if you did not see this model before

          # fit the model
          fit <- fit.function(item_assignment = item_assignment, 
                              item_names = item_names,
                              dat = data,
                              verbose = verbose,
                              fit_crit = fit_crit,
                              logistic_weights = logistic_weights,
                              nu_weights = nu_weights,
                              nu_min = nu_min,
                              debug_fit_mode = debug_fit_mode,
                              ...) #evaluate model
          
          
          # save items - factor allocation, seed, iteration, number of factors, "age" and quality in solutions
          n_nest_fac <- max(item_assignment)
          
                   
          new_sol <- c(unlist(item_assignment), 
                       seed = seed, 
                       iteration = iter, 
                       n_nest_fac = n_nest_fac, 
                       age = 0, 
                       fit)
        
        }
        new_sol
      }, simplify = "matrix")
      
    } else if (parallel){
      clusterExport(cl = myCL, varlist = ls(envir = environment()), envir = environment())
      clusterExport(cl = myCL, varlist = ls(envir = .GlobalEnv), envir = .GlobalEnv)
      
      new_solutions <- parSapply(myCL, 1:nrow(allBees), function(i_bee){
        
        if (verbose) print(allBees[i_bee,])
        
        item_assignment_orig <- allBees[i_bee, item_names]
        
        # change model depending on the role of the current bee
        if (allBees$role[i_bee] == "scout"){
          item_assignment <- scout.bee(item_assignment = item_assignment_orig,
                                       max_nest_fac = max_nest_fac)
          
          
        } else if (allBees$role[i_bee] == "onlooker") {
          
          item_assignment <- onlooker.bee(item_assignment = item_assignment_orig, 
                                          max_nest_fac = max_nest_fac)
          
        } else {stop("Fatal error: Bee swarm out of control.")} # if an unknown role occurs...
        
        #checks if this model has been tested before
        compare_item_assignment <- apply(solutions[, item_names], 1, function(i_solution) identical(item_assignment, i_solution))

        if (any(compare_item_assignment)){ # if you saw this item configuration before
          if (verbose) print ("Empty Flower")
          new_sol <- solutions[max(which(compare_item_assignment)),]
          new_sol$age <- new_sol$age + 1 
        } else { # if you did not see this model before
          
        
          fit <- fit.function(item_assignment = item_assignment, 
                              item_names = item_names,
                              dat = data,
                              verbose = verbose,
                              fit_crit = fit_crit,
                              logistic_weights = logistic_weights,
                              nu_weights = nu_weights,
                              nu_min = nu_min,
                              debug_fit_mode = debug_fit_mode,
                              ...
                              ) #evaluate model
          
          
          # save items - factor allocation, seed, iteration, number of factors, "age" and quality in solutions
          n_nest_fac <- max(item_assignment)
          
      
          new_sol <- c(unlist(item_assignment), 
                       seed = seed, 
                       iteration = iter, 
                       n_nest_fac = n_nest_fac, 
                       age = 0, 
                       fit)
    
          new_sol
        }
      }, simplify = "matrix")
    }
 
 
    # Transpose solutions for convenience purposes
    new_solutions <- data.frame(t(new_solutions))
    

    
    if (plot_nectar) {
      new_solutions2 <- data.frame(new_solutions) # create a work copy for plotting purposes
     # add necessary information for ggplot 
      new_solutions2$iter <- iter 
      new_solutions2$role <- as.factor(allBees$role)
      new_solutions2$n_nest_fac <- as.factor(new_solutions2$n_nest_fac)
      
      # plot and suppress NA-warnings (e.g., when fit is not converged)
      suppressWarnings(expr = {
        conv_plot <- conv_plot + 
          geom_jitter(data = new_solutions2, 
                      aes(x = iter,
                          y = fit_overall,
                          color = n_nest_fac,
                          shape = role),
                      width = plot_list$jitter_width,
                      alpha = plot_list$alpha,
                      size = plot_list$size) 
       if(!cluster_mode) plot(conv_plot)
      }) 
      
    } #plot quality of new solution
    

    # Sort solutions by nu so that the top_best are really on top
    
    if (all(is.na(new_solutions$fit_overall))) stop("All initial fit values are NA. Please check fit-function. Try setting verbose = TRUE or debug_fit_mode = TRUE.")


    new_solutions <- new_solutions[order(new_solutions$fit_overall, decreasing = TRUE),]
    
    if (max(new_solutions$fit_overall, na.rm = TRUE) > best_fit) { # if we have a new winner...
      best_fit <- max(new_solutions$fit_overall, na.rm = TRUE)
      best_solution <- new_solutions[1,]
      counter <- 0 # reset counter to zero because it will be counted +1 below
    } 
  
  
    # merge old and new solutions
    solutions <- rbind(solutions, new_solutions)
    solutions <- solutions[order(solutions$fit_overall, decreasing = T),] #sort by quality
    
    counter <- counter + 1 #increase count timer
    iter <- iter + 1   #increase overall iteration counter
    end <- Sys.time() 

  } #end search
  
  if (parallel) stopCluster(myCL)
  
  # write out final solutions
  write.table(solutions, summaryfile, sep=";",row.names = F, col.names = T, append = F)
  write.table(best_solution, summaryfile_fin, sep=";",row.names = F, col.names = T, append = F)
  
  # return conv_plot
  if(plot_nectar) {
    conv_plot + xlim(0,iter)
  } else {
      NULL
    }
  
}


############################################# End BSO #################################################

