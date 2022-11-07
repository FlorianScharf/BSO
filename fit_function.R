################## fit function #############################
#### This function contains model estimation and fit evaluation according 
#### to the custom criteria.
#### Please note that the fit function is highly customized for the bifactor model
#### investigated here and needs to be adapted fundamentally
#### in case BSO shall be applied to any other optimization problem.

fit.function = function(item_assignment, # a vector of item assignments as detailed in the section "Model Specification Search in Bifactor Modeling"
                        dat, # a data frame containing the data set
                        item_names, # a vector of item names (length must equal the length of item assignments)
                        fit_crit = c("cfi", "rmsea", "min_omega3", "min_loading"), # a vector of fit criteria to be extracted
                        # list (length equal to fit_crit) with weights for the logistic transformation (see section The Optimization Function)
                        logistic_weights = list(c(d = 0.9, a = 70),
                                                c(d = 0.06, a = -70),
                                                c(d = 0.4, a = 70),
                                                c(d = 0.33, a = 70)),
                        nu_weights = c(1,1,1,1), # a vector of weights for the aggregation function (see section The Optimization Function)
                        nu_min = 10^(-4), # threshold value, if the nu value of any fit criterion is below this threshold, nu
                        # will be returned as zero to prevent BSO from optimizing a subset of the fit criteria while neglecting
                        # one criterion completely
                        verbose = FALSE, # return information in console
                        return_fit = FALSE, # return fit object for debugging
                        debug_fit_mode = FALSE, # option to pass through error messages and warnings to help debugging the fit function
                        ... # arguments to be passed to lavaan
){
  
  
  
  tryCatch( { # not all model fits will converge, this prevents BSO from interrupting in those cases

    # show item assignments in verbose mode
    if (verbose) print(item_assignment)
    
    # the highest number among the item assignments represents the number of nested factors
    n_nest_fac = max(item_assignment)
    
    # write null model: explicitly model all variances so that the covariance matrix 
    # inside the lavaan object is complete irrespective of the chosen items
    null_model <- paste0(item_names,"~~",item_names, collapse = "\n")
    
    
    #### bifactor model definition ####
    # write model: modifier * general factor
    model <- paste0("g =~ ", 
                    paste0("lower(0) * ", 
                           item_names[item_assignment != -1], collapse = " + "))
    
    # write model: modifier * nested factor(s)
    for(fac in 1:n_nest_fac) {
      model <- paste0(model, "\n", "f", fac, " =~ ",
                      paste0("lower(0) * ",
                             item_names[item_assignment==fac], collapse = " + "))
      
    # Note that we restricted all factor loadings to be positive. 
    }
    #### Mpdel fitting ####

    # write current model in a text file in verbose mode
    if (verbose) writeLines(model, "Model_tmp.txt")
    
    # print current model to console in verbose mode
    if (verbose) print(model)
    
    # Fit Null Model
    fit_null <- cfa(model = null_model, 
                    data = dat, 
                    se = "none", # do not compute standard errors to speed up computation
                    ...) # additional arguments can be passed to lavaan from the main function (see lavOptions)
        

    
    # Fit current candidate model
    fit_model <- cfa(model = model, 
                     data = dat, 
                     std.lv = T, # standardize latent variables
                     se = "none", # do not compute standard errors to speed up computation
                     orthogonal = TRUE, # we assume orthogonal factors
                     ...) # additional arguments can be passed to lavaan from the main function (see lavOptions)
    
    
    #### Fit criterion evaluation ####
    ## In this section, we compute a number of fit criteria
    ## depending on the input to the fit criterion arguments
    
    # Preallocate vector for nectar values
    fit_inds <- c()
    
    # All fit indices in fit_crit which are computed in lavaan will be taken from lavaan
    fit_inds <- fitMeasures(fit_model, fit_crit)
    
    ### Compute reliability of nested factors
    rel_tmp <- semTools::compRelSEM(fit_model)
    rel <- rep(NA, max_nest_fac + 1) # length of the vector should be the same irrespective of the number of factors
    names(rel) <- c("g", paste0("f", 1:max_nest_fac))
    
    # sort rel values into results object
    rel[match(names(rel_tmp),names(rel))] <- rel_tmp
    ###
    
    ### Extract minimal reliability across all nested factors
    if ("min_omega3" %in% fit_crit){
      # take omega3 of all factors except the general factor
      # and extract the minimal value
      fit_inds <- c(fit_inds, min_omega3 = min(rel[-1], na.rm = TRUE))
    }
    ###
    
    ### Extract minimal (standardized) factor loading across all factors
    if ("min_loading" %in% fit_crit){
      # extract standardized loadings but only consider nested factors
      lambda_std <- lavInspect(fit_model, what = "std")$lambda[,-1]
      lambda_free <- (lavInspect(fit_model, what = "free")$lambda != 0)[,-1]
      
      # only consider estimated parameters
      lambda_std <- lambda_std[lambda_free]
      fit_inds <- c(fit_inds, min_loading = min(abs(lambda_std)))
    }
    ###
    
    ### Compute relative proportion of explained variance across all items
    if ("expl_var" %in% fit_crit){
      
      item_cov <-  lavInspect(fit_null, what = "sampstat")$cov
      
      item_resid_var <- diag(lavInspect(fit_model, what = "theta"))
      
      expl_var <- sum(diag(item_cov) - item_resid_var)/ sum(diag(item_cov))
      
      fit_inds <- c(fit_inds, expl_var = expl_var)
    }
    ###
    
    ### Extract relative proportion of included items
    if ("n_items" %in% fit_crit) {
      fit_inds <- c(fit_inds, n_items = sum(item_assignment > -1)/length(item_assignment))
    }
    ###

    
    
    #### Output preparation ####

    
    ### sort fit_inds according to order in argument for correct weighting
    fit_inds <- fit_inds[match(names(fit_inds), fit_crit)]
    
    # save detailed fit information
    fit_detailed <- round(c(fit_inds, 
                            rel)
                          , 3)
    
    # in verbose mode print detailed fit information to console
    if (verbose) print(fit_detailed)    
    
    
    # logistic transformation to bring all fit_inds on the same scale
    nu <- sapply(fit_crit, function(i_crit){
      index <- which(fit_crit == i_crit)
      out <- logistic(fit_inds[index], 
                      d = logistic_weights[[index]]["d"],
                      a = logistic_weights[[index]]["a"])
      names(out) <- NULL
      out
    })
    
    # aggregate nectar level nu according to vector of weights
    # set nu to zero if any value in nu is below nu_min to avoid 
    # neglecting of single fit criteria
    fit_overall <- ifelse(all(nu >= nu_min), yes = nu %*% nu_weights, no = 0)
    
    
    # if verbose print overall fit value
    if (verbose) print(fit_overall)
    
    # if return_fit (for debugging) return the lavaan fit object of the current model
    # otherwise return the fit values to be passed on to the best model table.
    if (return_fit) {
      fit_model 
    } else {
      c(fit_overall = fit_overall, nu = nu, fit_detailed)
    }
    
    ####
    
#### FAIL-SAFEs: Any warning or error will result in the fit function returning NA values
#### when the fit function is altered and a bug is introduced, all models may return NAs
#### in which case BSO returns a corresponding error message.     
  }, warning = function(w){
    if(debug_fit_mode){stop(w)}
    out_names <- c("fit_overall", paste0("nu.",fit_crit), fit_crit, c("g",paste0("f", 1:max_nest_fac)))
    out <- rep(NA,length(out_names))
    names(out) <- out_names
    out
  },
  error = function(e) {
    if(debug_fit_mode){stop(e)}
    out_names <- c("fit_overall", paste0("nu.",fit_crit), fit_crit, c("g",paste0("f", 1:max_nest_fac)))
    out <- rep(NA,length(out_names))
    names(out) <- out_names
    out
  })
  
  
} #end fit function
