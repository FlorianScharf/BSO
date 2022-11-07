################## BSO subfunctions #############################
#### This function contains a collection of sub-functions 
#### which are utilized during BSO.


#### Initilization function ####
# Generates the random start models in the first iteration
initialize.bees <- function(item_names, # a vector of item names (length must equal the length of item assignments)
                            data, # a data frame containing the data set
                            max_iter, # maximum number of iterations
                            n_bees, # number of bees
                            n_items = length(item_names), # number of items (inferred from name vector)
                            verbose = FALSE, # verbose mode (passed to fit function)
                            debug_fit_mode = FALSE, # debug mode (passed to fit function)
                            fit_crit = c("cfi", "rmsea", "min_omega3", "min_loading"), # see fit-function()
                            logistic_weights, # see fit-function()
                            nu_weights, # see fit-function()
                            nu_min, # see fit-function()
                            n_start_items = length(item_names), # How many items should be in the start model? (allows start with a subselection of items)
                            balance_n_fac = FALSE, # if TRUE, the start models are sampled such that the number of factors is uniformly distributed
                            ...){ # further arguments to be passed to lavaan
  
  
  # item_assignment is a vector with numbers with the following coding scheme:
  # -1                ... Item is not in the model
  # 0                 ... Item loads on the general factor only
  # 1:max_nest_fac    ... Item loads on this nested factor (and the general factor)
  
  
  
  if (!balance_n_fac){ 
    # if we start with a subset of the items, compute how many are excluded
    n_items_excluded <- n_items - n_start_items
    
    # return error for impossible argument combinations
    if(n_items_excluded < 0 | (n_items - n_items_excluded) < 0) stop("Invalid value for n_start_items")
    
    # assign random item allocations
    item_assignment <- rep(-1, n_items_excluded)
    item_assignment <- c(item_assignment, resamp(0:max_nest_fac, n_items - n_items_excluded, replace = TRUE)) #generate models
    
    # Note that across all start models this procedure will result in a higher number of models
    # with fewer factors (because these are more frequent among all possible random assignments)

  } else {
    
    # This approach tries to balance the number of factors across all start models. 
    
    # draw a random number of factors first
    n_fac_tmp <- sample(min_nest_fac:max_nest_fac, size = 1)
    item_assignment <- c(rep(1:n_fac_tmp, each = 3)) # place at least 3 items for each factor
    
    # the other items can be assigned to any of the factors 
    # items_in_model <- n_items - length(item_assignment)
    
    n_items_excluded <- n_items - n_start_items
      n_items_remaining <- n_items - length(item_assignment) - n_items_excluded
      # return error for impossible argument combinations
      if(n_items_excluded < 0 | n_items_remaining < 0) stop("Invalid value for n_start_items")
    item_assignment <- c(item_assignment, rep(-1, n_items_excluded), sample(0:n_fac_tmp, size = n_items_remaining, replace = TRUE))
    
    # randomize order
    item_assignment <- sample(item_assignment, size = length(item_assignment))
  }
  
 
  # remove nested factors with less than 3 items
  item_counts <- table(item_assignment[item_assignment %in% 1:max_nest_fac])
  fac_2b_removed <- as.numeric(names(item_counts)[item_counts < 3])
  # remove the factor by moving these items to the general factor only
  item_assignment[item_assignment %in% fac_2b_removed] <- 0
  
  # example: it is possible that a factor 5 occurs but factor 2 was removed
  # in that case, the factors should be relabeled 1:4
  # which items load on nested factors?
  item_assignment <- renew_factor_numbers(item_assignment = item_assignment)
  
  # Name assignment so that the names are exported
  names(item_assignment) <- item_names
  
  # if there are items loading on g only, the max fac is the general factor
  # otherwise it is 
  n_nest_fac <- max(item_assignment)
  
  # fit the the model
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
  
  # save items - factor allocation, seed, iteration, number of factors, "age" and quality in solution object
  new_sol <- c(item_assignment, 
               seed = seed, 
               iteration = 0, 
               n_nest_fac = n_nest_fac, 
               age = 0, 
               fit)
  new_sol 
}


#### Scout bee function ####
scout.bee <- function(item_assignment, max_nest_fac){
  n_items_per_fac <- table(factor(item_assignment, levels = 1:max_nest_fac))
  n_nest_fac   <- max(item_assignment)
  
  # Check which major changes are allowed and choose one randomly 
  allowed_operations <- check.operations.scouts(item_assignment = item_assignment,
                                                max_nest_fac = max_nest_fac)

  chosen_operation <- resamp(allowed_operations, size = 1) 

  switch (chosen_operation,
          { # add new factor
            free_items <- which(item_assignment <= 0) # Which items are available?
            n_new_items <- resamp(3:length(free_items), size = 1) # How many items should the new factor have?
            new_fac_items <- resamp(free_items, size = n_new_items) # Choose items randomly
            item_assignment[new_fac_items] <- n_nest_fac + 1}, # Assign these items to the new factor
          { #split factors
            big_fac <- which(n_items_per_fac >= 6) # Which factors could be split?
            split_fac <- resamp(big_fac, size = 1) # Decide which factor should be split
            split_fac_items <- which(item_assignment == split_fac) # Which items belong to the factors?
            n_new_fac_items <- resamp(3:(length(split_fac_items) - 3), size = 1) # Split the factor randomly such that each factor has at least 3 items
            new_fac_items <- resamp(split_fac_items, size = n_new_fac_items) # Sample items for new factor
            
            item_assignment[new_fac_items]<- n_nest_fac + 1},# Assign these items to the new factor
          { #remove factor
            old_fac <- resamp(1:n_nest_fac, size = 1)
            item_assignment[item_assignment == old_fac] <- 0 }, # Assign these items to the general factor
          { #merge factors
            old_fac <- resamp(1:n_nest_fac, size = 2)
            item_assignment[item_assignment == old_fac[1]] <- old_fac[2] }
  )
  
  # example: it is possible that a factor 5 occurs but factor 2 was removed
  # in that case, the factors should be relabeled 1:4
  # which items load on nested factors?
  item_assignment <- renew_factor_numbers(item_assignment = item_assignment)   
  item_assignment
}

#### Onlooker bee function ####
onlooker.bee <- function(item_assignment, max_nest_fac){
  n_items_per_fac <- table(factor(item_assignment,levels=1:max_nest_fac))
  n_nest_fac   <- max(item_assignment)

  # Check which major changes are allowed and choose one randomly
  allowed_operations <- check.operations.onlookers(item_assignment = item_assignment,
                                                max_nest_fac = max_nest_fac)

  chosen_operation <- resamp(allowed_operations, size = 1) 

  switch (chosen_operation,
          { # add item to nested factor
            which_item    <- resamp(which(item_assignment <= 0), size = 1) # look for items not assigned to a nested factor
            to_which_factor   <- resamp(0:n_nest_fac, 1) # sample a new factor for this item
            item_assignment[which_item] <- to_which_factor}, # change assignment
          { # remove item from nested factor
            which_factor   <- resamp(which(n_items_per_fac > 3), size = 1) # look for a factor with more than the minimal item count
            which_item   <- resamp(which(item_assignment == which_factor), size = 1) # sample any item from this factor
            item_assignment[which_item] <- 0}, # assign item to general factor
          { # swap item between nested factors
            item_1 <- resamp(which(item_assignment > 0), size = 1) # sample any item assigned to a nested factor
            item_2 <- resamp(which((item_assignment > 0) & (item_assignment != item_assignment[[item_1]])), size = 1) # sample any item from another nested factor
            item_assignment[c(item_1, item_2)] <- item_assignment[c(item_2, item_1)]},
          { # delete item from item pool
            candidate_factors   <- c(0, which(n_items_per_fac > 3)) # look for a factor with more than the minimal item count 
            which_item   <- resamp(which(item_assignment %in% candidate_factors), size = 1) # sample any item from this factor or pick an item from the general factor
            item_assignment[which_item] <- -1} # remove item from item pool completely
  )
  item_assignment # return
}

#### Scout bee operation selection ####
# decides which scout operations are possible based on an item assignment and 
# the number of factors
check.operations.scouts <- function(item_assignment, max_nest_fac){
  
  n_items_per_fac <- table(factor(item_assignment,levels=1:max_nest_fac))
  n_nest_fac   <- max(item_assignment)
  
  operations <- c(
    add_fac_allowed = (n_nest_fac < max_nest_fac) & (sum(item_assignment <= 0) >= 3),     #can scouts add factors?
    spl_fac_allowed = n_nest_fac < max_nest_fac & any(n_items_per_fac >= 6),  #can scouts split factors?
    rmv_fac_allowed = n_nest_fac > min_nest_fac,                      #can scouts remove factors?
    mer_fac_allowed = n_nest_fac > min_nest_fac                      #can scouts merge factors?
  )
   which(operations)
}
#### Onlooker bee operation selection ####
# decides which onlooker operations are possible based on an item assignment and 
# the number of factors
check.operations.onlookers <- function(item_assignment, max_nest_fac){
  
  n_items_per_fac <- table(factor(item_assignment,levels=1:max_nest_fac))
  n_nest_fac   <- max(item_assignment)
  
  operations <- c(
    add_item_allowed = any(item_assignment <= 0),         # can onlookers add item-factor allocations?
    rmv_item_allowed = any(n_items_per_fac > 3),   # can onlookers remove item-nested-factor allocations?
    swap_item_allowed = n_nest_fac > 1,         # can onlookers swap item-nested-factor allocations?
    del_item_allowed = any(item_assignment != -1) # can onlookers delete items from the model?
  )
  
   which(operations)
}


#### Minor helper functions ####
# ensure correct factor numbering
renew_factor_numbers <- function(item_assignment){
  # example: it is possible that a factor 5 occurs but factor 2 was removed
  # in that case, the factors should be relabeled 1:4
  # which items load on nested factors?
  items_nested <- item_assignment > 0
  item_assignment[items_nested] <- match(item_assignment[items_nested], sort(unique(item_assignment[items_nested])))
  item_assignment  
}

# sample function that will also sample correctly from vector length = 1
resamp <- function(y,...){if(length(y)==1) y else sample(y,...)} 

# debugging function to interrupt the script at any point
# and copy all temporary objects into the global environment
allglobal <- function(tmp_env = environment()) {
  lss <- ls(envir = tmp_env)
  for (i in lss) {
    assign(i, get(i, envir = tmp_env), envir = globalenv())
  }
}
