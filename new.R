library(stochtree)
bart_new <- function(X_train, y_train, W_train = NULL, group_ids_train = NULL, 
                 rfx_basis_train = NULL, X_test = NULL, W_test = NULL, 
                 group_ids_test = NULL, rfx_basis_test = NULL, 
                 num_gfr = 5, num_burnin = 0, num_mcmc = 100, 
                 params = list()) {
    # Extract BART parameters
    bart_params <- preprocessBartParams(params)
    cutpoint_grid_size <- bart_params$cutpoint_grid_size
    sigma_leaf_init <- bart_params$sigma_leaf_init
    alpha_mean <- bart_params$alpha_mean
    beta_mean <- bart_params$beta_mean
    min_samples_leaf_mean <- bart_params$min_samples_leaf_mean
    max_depth_mean <- bart_params$max_depth_mean
    alpha_variance <- bart_params$alpha_variance
    beta_variance <- bart_params$beta_variance
    min_samples_leaf_variance <- bart_params$min_samples_leaf_variance
    max_depth_variance <- bart_params$max_depth_variance
    a_global <- bart_params$a_global
    b_global <- bart_params$b_global
    a_leaf <- bart_params$a_leaf
    b_leaf <- bart_params$b_leaf
    a_forest <- bart_params$a_forest
    b_forest <- bart_params$b_forest
    variance_scale <- bart_params$variance_scale
    sigma2_init <- bart_params$sigma2_init
    variance_forest_init <- bart_params$variance_forest_init
    pct_var_sigma2_init <- bart_params$pct_var_sigma2_init
    pct_var_variance_forest_init <- bart_params$pct_var_variance_forest_init
    variable_weights_mean <- bart_params$variable_weights_mean
    variable_weights_variance <- bart_params$variable_weights_variance
    num_trees_mean <- bart_params$num_trees_mean
    num_trees_variance <- bart_params$num_trees_variance
    sample_sigma_global <- bart_params$sample_sigma_global
    sample_sigma_leaf <- bart_params$sample_sigma_leaf
    random_seed <- bart_params$random_seed
    keep_burnin <- bart_params$keep_burnin
    keep_gfr <- bart_params$keep_gfr
    verbose <- bart_params$verbose
    
    # Determine whether conditional mean, variance, or both will be modeled
    if (num_trees_variance > 0) include_variance_forest = T
    else include_variance_forest = F
    if (num_trees_mean > 0) include_mean_forest = T
    else include_mean_forest = F
    
    # Set the variance forest priors if not set
    if (include_variance_forest) {
        a_0 <- 1.5
        if (is.null(a_forest)) a_forest <- num_trees_variance / (a_0^2) + 0.5
        if (is.null(b_forest)) b_forest <- num_trees_variance / (a_0^2)
    } else {
        a_forest <- 1.
        b_forest <- 1.
    }
    
    # Override tau sampling if there is no mean forest
    if (!include_mean_forest) sample_sigma_leaf <- F
    
    # Variable weight preprocessing (and initialization if necessary)
    if (include_mean_forest) {
        if (is.null(variable_weights_mean)) {
            variable_weights_mean = rep(1/ncol(X_train), ncol(X_train))
        }
        if (any(variable_weights_mean < 0)) {
            stop("variable_weights_mean cannot have any negative weights")
        }
    }
    if (include_variance_forest) {
        if (is.null(variable_weights_variance)) {
            variable_weights_variance = rep(1/ncol(X_train), ncol(X_train))
        }
        if (any(variable_weights_variance < 0)) {
            stop("variable_weights_variance cannot have any negative weights")
        }
    }
    
    # Preprocess covariates
    if ((!is.data.frame(X_train)) && (!is.matrix(X_train))) {
        stop("X_train must be a matrix or dataframe")
    }
    if (!is.null(X_test)){
        if ((!is.data.frame(X_test)) && (!is.matrix(X_test))) {
            stop("X_test must be a matrix or dataframe")
        }
    }
    if ((ncol(X_train) != length(variable_weights_mean)) && (include_mean_forest)) {
        stop("length(variable_weights_mean) must equal ncol(X_train)")
    }
    if ((ncol(X_train) != length(variable_weights_variance)) && (include_variance_forest)) {
        stop("length(variable_weights_variance) must equal ncol(X_train)")
    }
    train_cov_preprocess_list <- preprocessTrainData(X_train)
    X_train_metadata <- train_cov_preprocess_list$metadata
    X_train <- train_cov_preprocess_list$data
    original_var_indices <- X_train_metadata$original_var_indices
    feature_types <- X_train_metadata$feature_types
    if (!is.null(X_test)) X_test <- preprocessPredictionData(X_test, X_train_metadata)
    
    # Update variable weights
    variable_weights_adj <- 1/sapply(original_var_indices, function(x) sum(original_var_indices == x))
    if (include_mean_forest) {
        variable_weights_mean <- variable_weights_mean[original_var_indices]*variable_weights_adj
    }
    if (include_variance_forest) {
        variable_weights_variance <- variable_weights_variance[original_var_indices]*variable_weights_adj
    }
    
    # Convert all input data to matrices if not already converted
    if ((is.null(dim(W_train))) && (!is.null(W_train))) {
        W_train <- as.matrix(W_train)
    }
    if ((is.null(dim(W_test))) && (!is.null(W_test))) {
        W_test <- as.matrix(W_test)
    }
    if ((is.null(dim(rfx_basis_train))) && (!is.null(rfx_basis_train))) {
        rfx_basis_train <- as.matrix(rfx_basis_train)
    }
    if ((is.null(dim(rfx_basis_test))) && (!is.null(rfx_basis_test))) {
        rfx_basis_test <- as.matrix(rfx_basis_test)
    }
    
    # Recode group IDs to integer vector (if passed as, for example, a vector of county names, etc...)
    has_rfx <- F
    has_rfx_test <- F
    if (!is.null(group_ids_train)) {
        group_ids_factor <- factor(group_ids_train)
        group_ids_train <- as.integer(group_ids_factor)
        has_rfx <- T
        if (!is.null(group_ids_test)) {
            group_ids_factor_test <- factor(group_ids_test, levels = levels(group_ids_factor))
            if (sum(is.na(group_ids_factor_test)) > 0) {
                stop("All random effect group labels provided in group_ids_test must be present in group_ids_train")
            }
            group_ids_test <- as.integer(group_ids_factor_test)
            has_rfx_test <- T
        }
    }
    
    # Data consistency checks
    if ((!is.null(X_test)) && (ncol(X_test) != ncol(X_train))) {
        stop("X_train and X_test must have the same number of columns")
    }
    if ((!is.null(W_test)) && (ncol(W_test) != ncol(W_train))) {
        stop("W_train and W_test must have the same number of columns")
    }
    if ((!is.null(W_train)) && (nrow(W_train) != nrow(X_train))) {
        stop("W_train and X_train must have the same number of rows")
    }
    if ((!is.null(W_test)) && (nrow(W_test) != nrow(X_test))) {
        stop("W_test and X_test must have the same number of rows")
    }
    if (nrow(X_train) != length(y_train)) {
        stop("X_train and y_train must have the same number of observations")
    }
    if ((!is.null(rfx_basis_test)) && (ncol(rfx_basis_test) != ncol(rfx_basis_train))) {
        stop("rfx_basis_train and rfx_basis_test must have the same number of columns")
    }
    if (!is.null(group_ids_train)) {
        if (!is.null(group_ids_test)) {
            if ((!is.null(rfx_basis_train)) && (is.null(rfx_basis_test))) {
                stop("rfx_basis_train is provided but rfx_basis_test is not provided")
            }
        }
    }
    
    # Fill in rfx basis as a vector of 1s (random intercept) if a basis not provided 
    has_basis_rfx <- F
    num_basis_rfx <- 0
    if (has_rfx) {
        if (is.null(rfx_basis_train)) {
            rfx_basis_train <- matrix(rep(1,nrow(X_train)), nrow = nrow(X_train), ncol = 1)
        } else {
            has_basis_rfx <- T
            num_basis_rfx <- ncol(rfx_basis_train)
        }
        num_rfx_groups <- length(unique(group_ids_train))
        num_rfx_components <- ncol(rfx_basis_train)
        if (num_rfx_groups == 1) warning("Only one group was provided for random effect sampling, so the 'redundant parameterization' is likely overkill")
    }
    if (has_rfx_test) {
        if (is.null(rfx_basis_test)) {
            if (!is.null(rfx_basis_train)) {
                stop("Random effects basis provided for training set, must also be provided for the test set")
            }
            rfx_basis_test <- matrix(rep(1,nrow(X_test)), nrow = nrow(X_test), ncol = 1)
        }
    }

    # Convert y_train to numeric vector if not already converted
    if (!is.null(dim(y_train))) {
        y_train <- as.matrix(y_train)
    }
    
    # Determine whether a basis vector is provided
    has_basis = !is.null(W_train)
    
    # Determine whether a test set is provided
    has_test = !is.null(X_test)

    # Standardize outcome separately for test and train
    y_bar_train <- mean(y_train)
    y_std_train <- sd(y_train)
    resid_train <- (y_train-y_bar_train)/y_std_train
    resid_train <- resid_train*sqrt(variance_scale)
    
    # Compute initial value of root nodes in mean forest
    init_val_mean <- mean(resid_train)

    # Calibrate priors for sigma^2 and tau
    if (is.null(sigma2_init)) sigma2_init <- pct_var_sigma2_init*var(resid_train)
    if (is.null(variance_forest_init)) variance_forest_init <- pct_var_variance_forest_init*var(resid_train)
    if (is.null(b_leaf)) b_leaf <- var(resid_train)/(2*num_trees_mean)
    if (has_basis) {
        if (ncol(W_train) > 1) {
            if (is.null(sigma_leaf_init)) sigma_leaf_init <- diag(var(resid_train)/(num_trees_mean), ncol(W_train))
            current_leaf_scale <- sigma_leaf_init
        } else {
            if (is.null(sigma_leaf_init)) sigma_leaf_init <- var(resid_train)/(num_trees_mean)
            current_leaf_scale <- as.matrix(sigma_leaf_init)
        }
    } else {
        if (is.null(sigma_leaf_init)) sigma_leaf_init <- var(resid_train)/(num_trees_mean)
        current_leaf_scale <- as.matrix(sigma_leaf_init)
    }
    current_sigma2 <- sigma2_init

    # Determine leaf model type
    if (!has_basis) leaf_model_mean_forest <- 0
    else if (ncol(W_train) == 1) leaf_model_mean_forest <- 1
    else if (ncol(W_train) > 1) leaf_model_mean_forest <- 2
    else stop("W_train passed must be a matrix with at least 1 column")

    # Set variance leaf model type (currently only one option)
    leaf_model_variance_forest <- 3
    
    # Unpack model type info
    if (leaf_model_mean_forest == 0) {
        output_dimension = 1
        is_leaf_constant = T
        leaf_regression = F
    } else if (leaf_model_mean_forest == 1) {
        stopifnot(has_basis)
        stopifnot(ncol(W_train) == 1)
        output_dimension = 1
        is_leaf_constant = F
        leaf_regression = T
    } else if (leaf_model_mean_forest == 2) {
        stopifnot(has_basis)
        stopifnot(ncol(W_train) > 1)
        output_dimension = ncol(W_train)
        is_leaf_constant = F
        leaf_regression = T
        if (sample_sigma_leaf) {
            stop("Sampling leaf scale not yet supported for multivariate leaf models")
        }
    }
    
    # Data
    if (leaf_regression) {
        forest_dataset_train <- createForestDataset(X_train, W_train)
        if (has_test) forest_dataset_test <- createForestDataset(X_test, W_test)
        requires_basis <- T
    } else {
        forest_dataset_train <- createForestDataset(X_train)
        if (has_test) forest_dataset_test <- createForestDataset(X_test)
        requires_basis <- F
    }
    outcome_train <- createOutcome(resid_train)
    
    # Random number generator (std::mt19937)
    if (is.null(random_seed)) random_seed = sample(1:10000,1,F)
    rng <- createRNG(random_seed)
    
    # Sampling data structures
    feature_types <- as.integer(feature_types)
    if (include_mean_forest) {
        current_forest_model_mean <- createForestModel(forest_dataset_train, feature_types, num_trees_mean, nrow(X_train), alpha_mean, beta_mean, min_samples_leaf_mean, max_depth_mean)
    }
    if (include_variance_forest) {
        forest_model_variance <- createForestModel(forest_dataset_train, feature_types, num_trees_variance, nrow(X_train), alpha_variance, beta_variance, min_samples_leaf_variance, max_depth_variance)
    }
    
    # Container of forest samples
    if (include_mean_forest) {
        current_forest_samples_mean <- createForestContainer(num_trees_mean, output_dimension, is_leaf_constant, FALSE)
    }
    if (include_variance_forest) {
        forest_samples_variance <- createForestContainer(num_trees_variance, 1, TRUE, TRUE)
    }
    
    # Random effects prior parameters
    if (has_rfx) {
        if (num_rfx_components == 1) {
            alpha_init <- c(1)
        } else if (num_rfx_components > 1) {
            alpha_init <- c(1,rep(0,num_rfx_components-1))
        } else {
            stop("There must be at least 1 random effect component")
        }
        xi_init <- matrix(rep(alpha_init, num_rfx_groups),num_rfx_components,num_rfx_groups)
        sigma_alpha_init <- diag(1,num_rfx_components,num_rfx_components)
        sigma_xi_init <- diag(1,num_rfx_components,num_rfx_components)
        sigma_xi_shape <- 1
        sigma_xi_scale <- 1
    }

    # Random effects data structure and storage container
    if (has_rfx) {
        rfx_dataset_train <- createRandomEffectsDataset(group_ids_train, rfx_basis_train)
        rfx_tracker_train <- createRandomEffectsTracker(group_ids_train)
        rfx_model <- createRandomEffectsModel(num_rfx_components, num_rfx_groups)
        rfx_model$set_working_parameter(alpha_init)
        rfx_model$set_group_parameters(xi_init)
        rfx_model$set_working_parameter_cov(sigma_alpha_init)
        rfx_model$set_group_parameter_cov(sigma_xi_init)
        rfx_model$set_variance_prior_shape(sigma_xi_shape)
        rfx_model$set_variance_prior_scale(sigma_xi_scale)
        rfx_samples <- createRandomEffectSamples(num_rfx_components, num_rfx_groups, rfx_tracker_train)
    }

    # Container of variance parameter samples
    num_samples <- num_gfr + num_burnin + num_mcmc
    if (sample_sigma_global) global_var_samples <- rep(0, num_samples)
    if (sample_sigma_leaf) leaf_scale_samples <- rep(0, num_samples)
    
    # Initialize the leaves of each tree in the mean forest
    if (include_mean_forest) {
        if (requires_basis) init_values_mean_forest <- rep(0., ncol(W_train))
        else init_values_mean_forest <- 0.
        current_forest_samples_mean$prepare_for_sampler(forest_dataset_train, outcome_train, current_forest_model_mean, leaf_model_mean_forest, init_values_mean_forest)
    }

    # Initialize the leaves of each tree in the variance forest
    if (include_variance_forest) {
        forest_samples_variance$prepare_for_sampler(forest_dataset_train, outcome_train, forest_model_variance, leaf_model_variance_forest, variance_forest_init)
    }

    num_tree_list <- c()
    num_tree_list<-c(num_tree_list,200)
    
    # Run GFR (warm start) if specified
    if (num_gfr > 0){
        gfr_indices = 1:num_gfr
        for (i in 1:num_gfr) {
            # Print progress
            if (verbose) {
                if ((i %% 10 == 0) || (i == num_gfr)) {
                    cat("Sampling", i, "out of", num_gfr, "XBART (grow-from-root) draws\n")
                }
            }
            
            if (include_mean_forest) {
                current_forest_model_mean$sample_one_iteration(
                    forest_dataset_train, outcome_train, current_forest_samples_mean, rng, feature_types, 
                    leaf_model_mean_forest, current_leaf_scale, variable_weights_mean, 
                    a_forest, b_forest, current_sigma2, cutpoint_grid_size, gfr = T, pre_initialized = T
                )
            }
            current_num_trees_mean <- rpois(1,200)
            current_forest_model_mean <- createForestModel(forest_dataset_train, feature_types, current_num_trees_mean, nrow(X_train), alpha_mean, beta_mean, min_samples_leaf_mean, max_depth_mean)
            current_forest_samples_mean <- createForestContainer(current_num_trees_mean, output_dimension, is_leaf_constant, FALSE)
            if (include_mean_forest) {
                if (requires_basis) init_values_mean_forest <- rep(0., ncol(W_train))
                else init_values_mean_forest <- 0.
                current_forest_samples_mean$prepare_for_sampler(forest_dataset_train, outcome_train, current_forest_model_mean, leaf_model_mean_forest, init_values_mean_forest)
            }
            num_tree_list <- c(num_tree_list, current_num_trees_mean)
            if (include_variance_forest) {
                forest_model_variance$sample_one_iteration(
                    forest_dataset_train, outcome_train, forest_samples_variance, rng, feature_types, 
                    leaf_model_variance_forest, current_leaf_scale, variable_weights_variance, 
                    a_forest, b_forest, current_sigma2, cutpoint_grid_size, gfr = T, pre_initialized = T
                )
            }
            if (sample_sigma_global) {
                global_var_samples[i] <- sample_sigma2_one_iteration(outcome_train, forest_dataset_train, rng, a_global, b_global)
                current_sigma2 <- global_var_samples[i]
            }
            if (sample_sigma_leaf) {
                leaf_scale_samples[i] <- sample_tau_one_iteration(current_forest_samples_mean, rng, a_leaf, b_leaf, i-1)
                current_leaf_scale <- as.matrix(leaf_scale_samples[i])
            }
            if (has_rfx) {
                rfx_model$sample_random_effect(rfx_dataset_train, outcome_train, rfx_tracker_train, rfx_samples, current_sigma2, rng)
            }
            if (include_mean_forest) {
                current_forest_model_mean$sample_one_iteration(
                    forest_dataset_train, outcome_train, current_forest_samples_mean, rng, feature_types, 
                    leaf_model_mean_forest, current_leaf_scale, variable_weights_mean, 
                    a_forest, b_forest, current_sigma2, cutpoint_grid_size, gfr = T, pre_initialized = T
                )
            }
        }
    }

    # Run MCMC
    if (num_burnin + num_mcmc > 0) {
        if (num_burnin > 0) {
            burnin_indices = (num_gfr+1):(num_gfr+num_burnin)
        }
        if (num_mcmc > 0) {
            mcmc_indices = (num_gfr+num_burnin+1):(num_gfr+num_burnin+num_mcmc)
        }
        for (i in (num_gfr+1):num_samples) {
            # Print progress
            if (verbose) {
                if (num_burnin > 0) {
                    if (((i - num_gfr) %% 100 == 0) || ((i - num_gfr) == num_burnin)) {
                        cat("Sampling", i - num_gfr, "out of", num_burnin, "BART burn-in draws\n")
                    }
                }
                if (num_mcmc > 0) {
                    if (((i - num_gfr - num_burnin) %% 100 == 0) || (i == num_samples)) {
                        cat("Sampling", i - num_burnin - num_gfr, "out of", num_mcmc, "BART MCMC draws\n")
                    }
                }
            }
            
            if (include_mean_forest) {
                current_forest_model_mean$sample_one_iteration(
                    forest_dataset_train, outcome_train, current_forest_samples_mean, rng, feature_types, 
                    leaf_model_mean_forest, current_leaf_scale, variable_weights_mean, 
                    a_forest, b_forest, current_sigma2, cutpoint_grid_size, gfr = F, pre_initialized = F
                )
            }
            current_num_trees_mean <- rpois(1,200)
            current_forest_model_mean <- createForestModel(forest_dataset_train, feature_types, current_num_trees_mean, nrow(X_train), alpha_mean, beta_mean, min_samples_leaf_mean, max_depth_mean)
            current_forest_samples_mean <- createForestContainer(current_num_trees_mean, output_dimension, is_leaf_constant, FALSE)
            if (include_mean_forest) {
                if (requires_basis) init_values_mean_forest <- rep(0., ncol(W_train))
                else init_values_mean_forest <- 0.
                current_forest_samples_mean$prepare_for_sampler(forest_dataset_train, outcome_train, current_forest_model_mean, leaf_model_mean_forest, init_values_mean_forest)
            }
            num_tree_list <- c(num_tree_list, current_num_trees_mean)
            if (include_variance_forest) {
                forest_model_variance$sample_one_iteration(
                    forest_dataset_train, outcome_train, forest_samples_variance, rng, feature_types, 
                    leaf_model_variance_forest, current_leaf_scale, variable_weights_variance, 
                    a_forest, b_forest, current_sigma2, cutpoint_grid_size, gfr = F, pre_initialized = T
                )
            }
            if (sample_sigma_global) {
                global_var_samples[i] <- sample_sigma2_one_iteration(outcome_train, forest_dataset_train, rng, a_global, b_global)
                current_sigma2 <- global_var_samples[i]
            }
            if (sample_sigma_leaf) {
                leaf_scale_samples[i] <- sample_tau_one_iteration(current_forest_samples_mean, rng, a_leaf, b_leaf, i-1)
                current_leaf_scale <- as.matrix(leaf_scale_samples[i])
            }
            if (has_rfx) {
                rfx_model$sample_random_effect(rfx_dataset_train, outcome_train, rfx_tracker_train, rfx_samples, current_sigma2, rng)
            }
            if (include_mean_forest) {
                current_forest_model_mean$sample_one_iteration(
                    forest_dataset_train, outcome_train, current_forest_samples_mean, rng, feature_types, 
                    leaf_model_mean_forest, current_leaf_scale, variable_weights_mean, 
                    a_forest, b_forest, current_sigma2, cutpoint_grid_size, gfr = F, pre_initialized = F
                )
            }
        }
    }
    
    if (include_mean_forest) {
        y_hat_train <- current_forest_samples_mean$predict(forest_dataset_train)*y_std_train/sqrt(variance_scale) + y_bar_train
        if (has_test) y_hat_test <- current_forest_samples_mean$predict(forest_dataset_test)*y_std_train/sqrt(variance_scale) + y_bar_train
    }
    
    # Variance forest predictions
    if (include_variance_forest) {
        sigma_x_hat_train <- forest_samples_variance$predict(forest_dataset_train)
        if (has_test) sigma_x_hat_test <- forest_samples_variance$predict(forest_dataset_test)
    }
    
    # Random effects predictions
    if (has_rfx) {
        rfx_preds_train <- rfx_samples$predict(group_ids_train, rfx_basis_train)*y_std_train/sqrt(variance_scale)
        y_hat_train <- y_hat_train + rfx_preds_train
    }
    if ((has_rfx_test) && (has_test)) {
        rfx_preds_test <- rfx_samples$predict(group_ids_test, rfx_basis_test)*y_std_train/sqrt(variance_scale)
        y_hat_test <- y_hat_test + rfx_preds_test
    }
    
    # Compute retention indices
    if (num_mcmc > 0) {
        keep_indices = mcmc_indices
        if (keep_gfr) keep_indices <- c(gfr_indices, keep_indices)
        if (keep_burnin) keep_indices <- c(burnin_indices, keep_indices)
    } else {
        if ((num_gfr > 0) && (num_burnin > 0)) {
            # Override keep_gfr = FALSE since there are no MCMC samples
            # Don't retain both GFR and burnin samples
            keep_indices = gfr_indices
        } else if ((num_gfr <= 0) && (num_burnin > 0)) {
            # Override keep_burnin = FALSE since there are no MCMC or GFR samples
            keep_indices = burnin_indices
        } else if ((num_gfr > 0) && (num_burnin <= 0)) {
            # Override keep_gfr = FALSE since there are no MCMC samples
            keep_indices = gfr_indices
        } else {
            stop("There are no samples to retain!")
        } 
    }
    keep_indices = 1
    # Subset forest and RFX predictions
    if (include_mean_forest) {
        y_hat_train <- y_hat_train[,keep_indices]
        if (has_test) y_hat_test <- y_hat_test[,keep_indices]
    }
    if (include_variance_forest) {
        sigma_x_hat_train <- sigma_x_hat_train[,keep_indices]
        if (has_test) sigma_x_hat_test <- sigma_x_hat_test[,keep_indices]
    }
    if (has_rfx) {
        rfx_preds_train <- rfx_preds_train[,keep_indices]
        if (has_test) rfx_preds_test <- rfx_preds_test[,keep_indices]
    }

    # Global error variance
    if (sample_sigma_global) sigma2_samples <- global_var_samples[keep_indices]*(y_std_train^2)/variance_scale
    
    # Leaf parameter variance
    if (sample_sigma_leaf) tau_samples <- leaf_scale_samples[keep_indices]
    
    # Rescale variance forest prediction by global sigma2 (sampled or constant)
    if (include_variance_forest) {
        if (sample_sigma_global) {
            sigma_x_hat_train <- sapply(1:length(keep_indices), function(i) sqrt(sigma_x_hat_train[,i]*sigma2_samples[i]))
            if (has_test) sigma_x_hat_test <- sapply(1:length(keep_indices), function(i) sqrt(sigma_x_hat_test[,i]*sigma2_samples[i]))
        } else {
            sigma_x_hat_train <- sqrt(sigma_x_hat_train*sigma2_init)*y_std_train/sqrt(variance_scale)
            if (has_test) sigma_x_hat_test <- sqrt(sigma_x_hat_test*sigma2_init)*y_std_train/sqrt(variance_scale)
        }
    }
    
    # Return results as a list
    # TODO: store variance_scale and propagate through predict function
    model_params <- list(
        "sigma2_init" = sigma2_init, 
        "sigma_leaf_init" = sigma_leaf_init,
        "a_global" = a_global,
        "b_global" = b_global, 
        "a_leaf" = a_leaf, 
        "b_leaf" = b_leaf,
        "a_forest" = a_forest, 
        "b_forest" = b_forest,
        "outcome_mean" = y_bar_train,
        "outcome_scale" = y_std_train, 
        "output_dimension" = output_dimension,
        "is_leaf_constant" = is_leaf_constant,
        "leaf_regression" = leaf_regression,
        "requires_basis" = requires_basis, 
        "num_covariates" = ncol(X_train), 
        "num_basis" = ifelse(is.null(W_train),0,ncol(W_train)), 
        "num_samples" = num_samples, 
        "num_gfr" = num_gfr, 
        "num_burnin" = num_burnin, 
        "num_mcmc" = num_mcmc, 
        "num_retained_samples" = length(keep_indices),
        "has_basis" = !is.null(W_train), 
        "has_rfx" = has_rfx, 
        "has_rfx_basis" = has_basis_rfx, 
        "num_rfx_basis" = num_basis_rfx, 
        "sample_sigma_global" = sample_sigma_global,
        "sample_sigma_leaf" = sample_sigma_leaf,
        "include_mean_forest" = include_mean_forest,
        "include_variance_forest" = include_variance_forest,
        "variance_scale" = variance_scale,
        "num_tree_list" = num_tree_list
    )
    result <- list(
        "model_params" = model_params, 
        "train_set_metadata" = X_train_metadata,
        "keep_indices" = keep_indices
    )
    if (include_mean_forest) {
        result[["mean_forests"]] = current_forest_samples_mean
        result[["y_hat_train"]] = y_hat_train
        if (has_test) result[["y_hat_test"]] = y_hat_test
    }
    if (include_variance_forest) {
        result[["variance_forests"]] = forest_samples_variance
        result[["sigma_x_hat_train"]] = sigma_x_hat_train
        if (has_test) result[["sigma_x_hat_test"]] = sigma_x_hat_test
    }
    if (sample_sigma_global) result[["sigma2_global_samples"]] = sigma2_samples
    if (sample_sigma_leaf) result[["sigma2_leaf_samples"]] = tau_samples
    if (has_rfx) {
        result[["rfx_samples"]] = rfx_samples
        result[["rfx_preds_train"]] = rfx_preds_train
        result[["rfx_unique_group_ids"]] = levels(group_ids_factor)
    }
    if ((has_rfx_test) && (has_test)) result[["rfx_preds_test"]] = rfx_preds_test
    class(result) <- "bartmodel"
    
    # Clean up classes with external pointers to C++ data structures
    if (include_mean_forest) rm(current_forest_model_mean)
    if (include_variance_forest) rm(forest_model_variance)
    rm(forest_dataset_train)
    if (has_test) rm(forest_dataset_test)
    if (has_rfx) rm(rfx_dataset_train, rfx_tracker_train, rfx_model)
    rm(outcome_train)
    rm(rng)
    
    return(result)
}

bart_test <- function(X_train, y_train, W_train = NULL, group_ids_train = NULL, 
                 rfx_basis_train = NULL, X_test = NULL, W_test = NULL, 
                 group_ids_test = NULL, rfx_basis_test = NULL, 
                 num_gfr = 5, num_burnin = 0, num_mcmc = 100, 
                 params = list()) {
    # Extract BART parameters
    bart_params <- preprocessBartParams(params)
    cutpoint_grid_size <- bart_params$cutpoint_grid_size
    sigma_leaf_init <- bart_params$sigma_leaf_init
    alpha_mean <- bart_params$alpha_mean
    beta_mean <- bart_params$beta_mean
    min_samples_leaf_mean <- bart_params$min_samples_leaf_mean
    max_depth_mean <- bart_params$max_depth_mean
    alpha_variance <- bart_params$alpha_variance
    beta_variance <- bart_params$beta_variance
    min_samples_leaf_variance <- bart_params$min_samples_leaf_variance
    max_depth_variance <- bart_params$max_depth_variance
    a_global <- bart_params$a_global
    b_global <- bart_params$b_global
    a_leaf <- bart_params$a_leaf
    b_leaf <- bart_params$b_leaf
    a_forest <- bart_params$a_forest
    b_forest <- bart_params$b_forest
    variance_scale <- bart_params$variance_scale
    sigma2_init <- bart_params$sigma2_init
    variance_forest_init <- bart_params$variance_forest_init
    pct_var_sigma2_init <- bart_params$pct_var_sigma2_init
    pct_var_variance_forest_init <- bart_params$pct_var_variance_forest_init
    variable_weights_mean <- bart_params$variable_weights_mean
    variable_weights_variance <- bart_params$variable_weights_variance
    num_trees_mean <- bart_params$num_trees_mean
    num_trees_variance <- bart_params$num_trees_variance
    sample_sigma_global <- bart_params$sample_sigma_global
    sample_sigma_leaf <- bart_params$sample_sigma_leaf
    random_seed <- bart_params$random_seed
    keep_burnin <- bart_params$keep_burnin
    keep_gfr <- bart_params$keep_gfr
    verbose <- bart_params$verbose
    
    # Determine whether conditional mean, variance, or both will be modeled
    if (num_trees_variance > 0) include_variance_forest = T
    else include_variance_forest = F
    if (num_trees_mean > 0) include_mean_forest = T
    else include_mean_forest = F
    
    # Set the variance forest priors if not set
    if (include_variance_forest) {
        a_0 <- 1.5
        if (is.null(a_forest)) a_forest <- num_trees_variance / (a_0^2) + 0.5
        if (is.null(b_forest)) b_forest <- num_trees_variance / (a_0^2)
    } else {
        a_forest <- 1.
        b_forest <- 1.
    }
    
    # Override tau sampling if there is no mean forest
    if (!include_mean_forest) sample_sigma_leaf <- F
    
    # Variable weight preprocessing (and initialization if necessary)
    if (include_mean_forest) {
        if (is.null(variable_weights_mean)) {
            variable_weights_mean = rep(1/ncol(X_train), ncol(X_train))
        }
        if (any(variable_weights_mean < 0)) {
            stop("variable_weights_mean cannot have any negative weights")
        }
    }
    if (include_variance_forest) {
        if (is.null(variable_weights_variance)) {
            variable_weights_variance = rep(1/ncol(X_train), ncol(X_train))
        }
        if (any(variable_weights_variance < 0)) {
            stop("variable_weights_variance cannot have any negative weights")
        }
    }
    
    # Preprocess covariates
    if ((!is.data.frame(X_train)) && (!is.matrix(X_train))) {
        stop("X_train must be a matrix or dataframe")
    }
    if (!is.null(X_test)){
        if ((!is.data.frame(X_test)) && (!is.matrix(X_test))) {
            stop("X_test must be a matrix or dataframe")
        }
    }
    if ((ncol(X_train) != length(variable_weights_mean)) && (include_mean_forest)) {
        stop("length(variable_weights_mean) must equal ncol(X_train)")
    }
    if ((ncol(X_train) != length(variable_weights_variance)) && (include_variance_forest)) {
        stop("length(variable_weights_variance) must equal ncol(X_train)")
    }
    train_cov_preprocess_list <- preprocessTrainData(X_train)
    X_train_metadata <- train_cov_preprocess_list$metadata
    X_train <- train_cov_preprocess_list$data
    original_var_indices <- X_train_metadata$original_var_indices
    feature_types <- X_train_metadata$feature_types
    if (!is.null(X_test)) X_test <- preprocessPredictionData(X_test, X_train_metadata)
    
    # Update variable weights
    variable_weights_adj <- 1/sapply(original_var_indices, function(x) sum(original_var_indices == x))
    if (include_mean_forest) {
        variable_weights_mean <- variable_weights_mean[original_var_indices]*variable_weights_adj
    }
    if (include_variance_forest) {
        variable_weights_variance <- variable_weights_variance[original_var_indices]*variable_weights_adj
    }
    
    # Convert all input data to matrices if not already converted
    if ((is.null(dim(W_train))) && (!is.null(W_train))) {
        W_train <- as.matrix(W_train)
    }
    if ((is.null(dim(W_test))) && (!is.null(W_test))) {
        W_test <- as.matrix(W_test)
    }
    if ((is.null(dim(rfx_basis_train))) && (!is.null(rfx_basis_train))) {
        rfx_basis_train <- as.matrix(rfx_basis_train)
    }
    if ((is.null(dim(rfx_basis_test))) && (!is.null(rfx_basis_test))) {
        rfx_basis_test <- as.matrix(rfx_basis_test)
    }
    
    # Recode group IDs to integer vector (if passed as, for example, a vector of county names, etc...)
    has_rfx <- F
    has_rfx_test <- F
    if (!is.null(group_ids_train)) {
        group_ids_factor <- factor(group_ids_train)
        group_ids_train <- as.integer(group_ids_factor)
        has_rfx <- T
        if (!is.null(group_ids_test)) {
            group_ids_factor_test <- factor(group_ids_test, levels = levels(group_ids_factor))
            if (sum(is.na(group_ids_factor_test)) > 0) {
                stop("All random effect group labels provided in group_ids_test must be present in group_ids_train")
            }
            group_ids_test <- as.integer(group_ids_factor_test)
            has_rfx_test <- T
        }
    }
    
    # Data consistency checks
    if ((!is.null(X_test)) && (ncol(X_test) != ncol(X_train))) {
        stop("X_train and X_test must have the same number of columns")
    }
    if ((!is.null(W_test)) && (ncol(W_test) != ncol(W_train))) {
        stop("W_train and W_test must have the same number of columns")
    }
    if ((!is.null(W_train)) && (nrow(W_train) != nrow(X_train))) {
        stop("W_train and X_train must have the same number of rows")
    }
    if ((!is.null(W_test)) && (nrow(W_test) != nrow(X_test))) {
        stop("W_test and X_test must have the same number of rows")
    }
    if (nrow(X_train) != length(y_train)) {
        stop("X_train and y_train must have the same number of observations")
    }
    if ((!is.null(rfx_basis_test)) && (ncol(rfx_basis_test) != ncol(rfx_basis_train))) {
        stop("rfx_basis_train and rfx_basis_test must have the same number of columns")
    }
    if (!is.null(group_ids_train)) {
        if (!is.null(group_ids_test)) {
            if ((!is.null(rfx_basis_train)) && (is.null(rfx_basis_test))) {
                stop("rfx_basis_train is provided but rfx_basis_test is not provided")
            }
        }
    }
    
    # Fill in rfx basis as a vector of 1s (random intercept) if a basis not provided 
    has_basis_rfx <- F
    num_basis_rfx <- 0
    if (has_rfx) {
        if (is.null(rfx_basis_train)) {
            rfx_basis_train <- matrix(rep(1,nrow(X_train)), nrow = nrow(X_train), ncol = 1)
        } else {
            has_basis_rfx <- T
            num_basis_rfx <- ncol(rfx_basis_train)
        }
        num_rfx_groups <- length(unique(group_ids_train))
        num_rfx_components <- ncol(rfx_basis_train)
        if (num_rfx_groups == 1) warning("Only one group was provided for random effect sampling, so the 'redundant parameterization' is likely overkill")
    }
    if (has_rfx_test) {
        if (is.null(rfx_basis_test)) {
            if (!is.null(rfx_basis_train)) {
                stop("Random effects basis provided for training set, must also be provided for the test set")
            }
            rfx_basis_test <- matrix(rep(1,nrow(X_test)), nrow = nrow(X_test), ncol = 1)
        }
    }

    # Convert y_train to numeric vector if not already converted
    if (!is.null(dim(y_train))) {
        y_train <- as.matrix(y_train)
    }
    
    # Determine whether a basis vector is provided
    has_basis = !is.null(W_train)
    
    # Determine whether a test set is provided
    has_test = !is.null(X_test)

    # Standardize outcome separately for test and train
    y_bar_train <- mean(y_train)
    y_std_train <- sd(y_train)
    resid_train <- (y_train-y_bar_train)/y_std_train
    resid_train <- resid_train*sqrt(variance_scale)
    
    # Compute initial value of root nodes in mean forest
    init_val_mean <- mean(resid_train)

    # Calibrate priors for sigma^2 and tau
    if (is.null(sigma2_init)) sigma2_init <- pct_var_sigma2_init*var(resid_train)
    if (is.null(variance_forest_init)) variance_forest_init <- pct_var_variance_forest_init*var(resid_train)
    if (is.null(b_leaf)) b_leaf <- var(resid_train)/(2*num_trees_mean)
    if (has_basis) {
        if (ncol(W_train) > 1) {
            if (is.null(sigma_leaf_init)) sigma_leaf_init <- diag(var(resid_train)/(num_trees_mean), ncol(W_train))
            current_leaf_scale <- sigma_leaf_init
        } else {
            if (is.null(sigma_leaf_init)) sigma_leaf_init <- var(resid_train)/(num_trees_mean)
            current_leaf_scale <- as.matrix(sigma_leaf_init)
        }
    } else {
        if (is.null(sigma_leaf_init)) sigma_leaf_init <- var(resid_train)/(num_trees_mean)
        current_leaf_scale <- as.matrix(sigma_leaf_init)
    }
    current_sigma2 <- sigma2_init

    # Determine leaf model type
    if (!has_basis) leaf_model_mean_forest <- 0
    else if (ncol(W_train) == 1) leaf_model_mean_forest <- 1
    else if (ncol(W_train) > 1) leaf_model_mean_forest <- 2
    else stop("W_train passed must be a matrix with at least 1 column")

    # Set variance leaf model type (currently only one option)
    leaf_model_variance_forest <- 3
    
    # Unpack model type info
    if (leaf_model_mean_forest == 0) {
        output_dimension = 1
        is_leaf_constant = T
        leaf_regression = F
    } else if (leaf_model_mean_forest == 1) {
        stopifnot(has_basis)
        stopifnot(ncol(W_train) == 1)
        output_dimension = 1
        is_leaf_constant = F
        leaf_regression = T
    } else if (leaf_model_mean_forest == 2) {
        stopifnot(has_basis)
        stopifnot(ncol(W_train) > 1)
        output_dimension = ncol(W_train)
        is_leaf_constant = F
        leaf_regression = T
        if (sample_sigma_leaf) {
            stop("Sampling leaf scale not yet supported for multivariate leaf models")
        }
    }
    
    # Data
    if (leaf_regression) {
        forest_dataset_train <- createForestDataset(X_train, W_train)
        if (has_test) forest_dataset_test <- createForestDataset(X_test, W_test)
        requires_basis <- T
    } else {
        forest_dataset_train <- createForestDataset(X_train)
        if (has_test) forest_dataset_test <- createForestDataset(X_test)
        requires_basis <- F
    }
    outcome_train <- createOutcome(resid_train)
    
    # Random number generator (std::mt19937)
    if (is.null(random_seed)) random_seed = sample(1:10000,1,F)
    rng <- createRNG(random_seed)
    
    # Sampling data structures
    feature_types <- as.integer(feature_types)
    if (include_mean_forest) {
        forest_model_mean <- createForestModel(forest_dataset_train, feature_types, num_trees_mean, nrow(X_train), alpha_mean, beta_mean, min_samples_leaf_mean, max_depth_mean)
    }
    if (include_variance_forest) {
        forest_model_variance <- createForestModel(forest_dataset_train, feature_types, num_trees_variance, nrow(X_train), alpha_variance, beta_variance, min_samples_leaf_variance, max_depth_variance)
    }
    
    # Container of forest samples
    if (include_mean_forest) {
        forest_samples_mean <- createForestContainer(num_trees_mean, output_dimension, is_leaf_constant, FALSE)
    }
    if (include_variance_forest) {
        forest_samples_variance <- createForestContainer(num_trees_variance, 1, TRUE, TRUE)
    }
    
    # Random effects prior parameters
    if (has_rfx) {
        if (num_rfx_components == 1) {
            alpha_init <- c(1)
        } else if (num_rfx_components > 1) {
            alpha_init <- c(1,rep(0,num_rfx_components-1))
        } else {
            stop("There must be at least 1 random effect component")
        }
        xi_init <- matrix(rep(alpha_init, num_rfx_groups),num_rfx_components,num_rfx_groups)
        sigma_alpha_init <- diag(1,num_rfx_components,num_rfx_components)
        sigma_xi_init <- diag(1,num_rfx_components,num_rfx_components)
        sigma_xi_shape <- 1
        sigma_xi_scale <- 1
    }

    # Random effects data structure and storage container
    if (has_rfx) {
        rfx_dataset_train <- createRandomEffectsDataset(group_ids_train, rfx_basis_train)
        rfx_tracker_train <- createRandomEffectsTracker(group_ids_train)
        rfx_model <- createRandomEffectsModel(num_rfx_components, num_rfx_groups)
        rfx_model$set_working_parameter(alpha_init)
        rfx_model$set_group_parameters(xi_init)
        rfx_model$set_working_parameter_cov(sigma_alpha_init)
        rfx_model$set_group_parameter_cov(sigma_xi_init)
        rfx_model$set_variance_prior_shape(sigma_xi_shape)
        rfx_model$set_variance_prior_scale(sigma_xi_scale)
        rfx_samples <- createRandomEffectSamples(num_rfx_components, num_rfx_groups, rfx_tracker_train)
    }

    # Container of variance parameter samples
    num_samples <- num_gfr + num_burnin + num_mcmc
    if (sample_sigma_global) global_var_samples <- rep(0, num_samples)
    if (sample_sigma_leaf) leaf_scale_samples <- rep(0, num_samples)
    
    # Initialize the leaves of each tree in the mean forest
    if (include_mean_forest) {
        if (requires_basis) init_values_mean_forest <- rep(0., ncol(W_train))
        else init_values_mean_forest <- 0.
        forest_samples_mean$prepare_for_sampler(forest_dataset_train, outcome_train, forest_model_mean, leaf_model_mean_forest, init_values_mean_forest)
    }

    # Initialize the leaves of each tree in the variance forest
    if (include_variance_forest) {
        forest_samples_variance$prepare_for_sampler(forest_dataset_train, outcome_train, forest_model_variance, leaf_model_variance_forest, variance_forest_init)
    }
    
    # Run GFR (warm start) if specified
    if (num_gfr > 0){
        gfr_indices = 1:num_gfr
        for (i in 1:num_gfr) {
            # Print progress
            if (verbose) {
                if ((i %% 10 == 0) || (i == num_gfr)) {
                    cat("Sampling", i, "out of", num_gfr, "XBART (grow-from-root) draws\n")
                }
            }
            
            if (include_mean_forest) {
                forest_model_mean$sample_one_iteration(
                    forest_dataset_train, outcome_train, forest_samples_mean, rng, feature_types, 
                    leaf_model_mean_forest, current_leaf_scale, variable_weights_mean, 
                    a_forest, b_forest, current_sigma2, cutpoint_grid_size, gfr = T, pre_initialized = T
                )
            }
            if (include_variance_forest) {
                forest_model_variance$sample_one_iteration(
                    forest_dataset_train, outcome_train, forest_samples_variance, rng, feature_types, 
                    leaf_model_variance_forest, current_leaf_scale, variable_weights_variance, 
                    a_forest, b_forest, current_sigma2, cutpoint_grid_size, gfr = T, pre_initialized = T
                )
            }
            if (sample_sigma_global) {
                global_var_samples[i] <- sample_sigma2_one_iteration(outcome_train, forest_dataset_train, rng, a_global, b_global)
                current_sigma2 <- global_var_samples[i]
            }
            if (sample_sigma_leaf) {
                leaf_scale_samples[i] <- sample_tau_one_iteration(forest_samples_mean, rng, a_leaf, b_leaf, i-1)
                current_leaf_scale <- as.matrix(leaf_scale_samples[i])
            }
            if (has_rfx) {
                rfx_model$sample_random_effect(rfx_dataset_train, outcome_train, rfx_tracker_train, rfx_samples, current_sigma2, rng)
            }
        }
    }
    
    # Run MCMC
    if (num_burnin + num_mcmc > 0) {
        if (num_burnin > 0) {
            burnin_indices = (num_gfr+1):(num_gfr+num_burnin)
        }
        if (num_mcmc > 0) {
            mcmc_indices = (num_gfr+num_burnin+1):(num_gfr+num_burnin+num_mcmc)
        }
        for (i in (num_gfr+1):num_samples) {
            # Print progress
            if (verbose) {
                if (num_burnin > 0) {
                    if (((i - num_gfr) %% 100 == 0) || ((i - num_gfr) == num_burnin)) {
                        cat("Sampling", i - num_gfr, "out of", num_burnin, "BART burn-in draws\n")
                    }
                }
                if (num_mcmc > 0) {
                    if (((i - num_gfr - num_burnin) %% 100 == 0) || (i == num_samples)) {
                        cat("Sampling", i - num_burnin - num_gfr, "out of", num_mcmc, "BART MCMC draws\n")
                    }
                }
            }
            
            if (include_mean_forest) {
                forest_model_mean$sample_one_iteration(
                    forest_dataset_train, outcome_train, forest_samples_mean, rng, feature_types, 
                    leaf_model_mean_forest, current_leaf_scale, variable_weights_mean, 
                    a_forest, b_forest, current_sigma2, cutpoint_grid_size, gfr = F, pre_initialized = T
                )
            }
            if (include_variance_forest) {
                forest_model_variance$sample_one_iteration(
                    forest_dataset_train, outcome_train, forest_samples_variance, rng, feature_types, 
                    leaf_model_variance_forest, current_leaf_scale, variable_weights_variance, 
                    a_forest, b_forest, current_sigma2, cutpoint_grid_size, gfr = F, pre_initialized = T
                )
            }
            if (sample_sigma_global) {
                global_var_samples[i] <- sample_sigma2_one_iteration(outcome_train, forest_dataset_train, rng, a_global, b_global)
                current_sigma2 <- global_var_samples[i]
            }
            if (sample_sigma_leaf) {
                leaf_scale_samples[i] <- sample_tau_one_iteration(forest_samples_mean, rng, a_leaf, b_leaf, i-1)
                current_leaf_scale <- as.matrix(leaf_scale_samples[i])
            }
            if (has_rfx) {
                rfx_model$sample_random_effect(rfx_dataset_train, outcome_train, rfx_tracker_train, rfx_samples, current_sigma2, rng)
            }
        }
    }
    
    # Mean forest predictions
    if (include_mean_forest) {
        y_hat_train <- forest_samples_mean$predict(forest_dataset_train)*y_std_train/sqrt(variance_scale) + y_bar_train
        if (has_test) y_hat_test <- forest_samples_mean$predict(forest_dataset_test)*y_std_train/sqrt(variance_scale) + y_bar_train
    }
    print(y_hat_train)
    # Variance forest predictions
    if (include_variance_forest) {
        sigma_x_hat_train <- forest_samples_variance$predict(forest_dataset_train)
        if (has_test) sigma_x_hat_test <- forest_samples_variance$predict(forest_dataset_test)
    }
    
    # Random effects predictions
    if (has_rfx) {
        rfx_preds_train <- rfx_samples$predict(group_ids_train, rfx_basis_train)*y_std_train/sqrt(variance_scale)
        y_hat_train <- y_hat_train + rfx_preds_train
    }
    if ((has_rfx_test) && (has_test)) {
        rfx_preds_test <- rfx_samples$predict(group_ids_test, rfx_basis_test)*y_std_train/sqrt(variance_scale)
        y_hat_test <- y_hat_test + rfx_preds_test
    }
    
    # Compute retention indices
    if (num_mcmc > 0) {
        keep_indices = mcmc_indices
        if (keep_gfr) keep_indices <- c(gfr_indices, keep_indices)
        if (keep_burnin) keep_indices <- c(burnin_indices, keep_indices)
    } else {
        if ((num_gfr > 0) && (num_burnin > 0)) {
            # Override keep_gfr = FALSE since there are no MCMC samples
            # Don't retain both GFR and burnin samples
            keep_indices = gfr_indices
        } else if ((num_gfr <= 0) && (num_burnin > 0)) {
            # Override keep_burnin = FALSE since there are no MCMC or GFR samples
            keep_indices = burnin_indices
        } else if ((num_gfr > 0) && (num_burnin <= 0)) {
            # Override keep_gfr = FALSE since there are no MCMC samples
            keep_indices = gfr_indices
        } else {
            stop("There are no samples to retain!")
        } 
    }
    
    # Subset forest and RFX predictions
    if (include_mean_forest) {
        y_hat_train <- y_hat_train[,keep_indices]
        if (has_test) y_hat_test <- y_hat_test[,keep_indices]
    }
    if (include_variance_forest) {
        sigma_x_hat_train <- sigma_x_hat_train[,keep_indices]
        if (has_test) sigma_x_hat_test <- sigma_x_hat_test[,keep_indices]
    }
    if (has_rfx) {
        rfx_preds_train <- rfx_preds_train[,keep_indices]
        if (has_test) rfx_preds_test <- rfx_preds_test[,keep_indices]
    }

    # Global error variance
    if (sample_sigma_global) sigma2_samples <- global_var_samples[keep_indices]*(y_std_train^2)/variance_scale
    
    # Leaf parameter variance
    if (sample_sigma_leaf) tau_samples <- leaf_scale_samples[keep_indices]
    
    # Rescale variance forest prediction by global sigma2 (sampled or constant)
    if (include_variance_forest) {
        if (sample_sigma_global) {
            sigma_x_hat_train <- sapply(1:length(keep_indices), function(i) sqrt(sigma_x_hat_train[,i]*sigma2_samples[i]))
            if (has_test) sigma_x_hat_test <- sapply(1:length(keep_indices), function(i) sqrt(sigma_x_hat_test[,i]*sigma2_samples[i]))
        } else {
            sigma_x_hat_train <- sqrt(sigma_x_hat_train*sigma2_init)*y_std_train/sqrt(variance_scale)
            if (has_test) sigma_x_hat_test <- sqrt(sigma_x_hat_test*sigma2_init)*y_std_train/sqrt(variance_scale)
        }
    }
    
    # Return results as a list
    # TODO: store variance_scale and propagate through predict function
    model_params <- list(
        "sigma2_init" = sigma2_init, 
        "sigma_leaf_init" = sigma_leaf_init,
        "a_global" = a_global,
        "b_global" = b_global, 
        "a_leaf" = a_leaf, 
        "b_leaf" = b_leaf,
        "a_forest" = a_forest, 
        "b_forest" = b_forest,
        "outcome_mean" = y_bar_train,
        "outcome_scale" = y_std_train, 
        "output_dimension" = output_dimension,
        "is_leaf_constant" = is_leaf_constant,
        "leaf_regression" = leaf_regression,
        "requires_basis" = requires_basis, 
        "num_covariates" = ncol(X_train), 
        "num_basis" = ifelse(is.null(W_train),0,ncol(W_train)), 
        "num_samples" = num_samples, 
        "num_gfr" = num_gfr, 
        "num_burnin" = num_burnin, 
        "num_mcmc" = num_mcmc, 
        "num_retained_samples" = length(keep_indices),
        "has_basis" = !is.null(W_train), 
        "has_rfx" = has_rfx, 
        "has_rfx_basis" = has_basis_rfx, 
        "num_rfx_basis" = num_basis_rfx, 
        "sample_sigma_global" = sample_sigma_global,
        "sample_sigma_leaf" = sample_sigma_leaf,
        "include_mean_forest" = include_mean_forest,
        "include_variance_forest" = include_variance_forest,
        "variance_scale" = variance_scale
    )
    result <- list(
        "model_params" = model_params, 
        "train_set_metadata" = X_train_metadata,
        "keep_indices" = keep_indices
    )
    if (include_mean_forest) {
        result[["mean_forests"]] = forest_samples_mean
        result[["y_hat_train"]] = y_hat_train
        if (has_test) result[["y_hat_test"]] = y_hat_test
    }
    if (include_variance_forest) {
        result[["variance_forests"]] = forest_samples_variance
        result[["sigma_x_hat_train"]] = sigma_x_hat_train
        if (has_test) result[["sigma_x_hat_test"]] = sigma_x_hat_test
    }
    if (sample_sigma_global) result[["sigma2_global_samples"]] = sigma2_samples
    if (sample_sigma_leaf) result[["sigma2_leaf_samples"]] = tau_samples
    if (has_rfx) {
        result[["rfx_samples"]] = rfx_samples
        result[["rfx_preds_train"]] = rfx_preds_train
        result[["rfx_unique_group_ids"]] = levels(group_ids_factor)
    }
    if ((has_rfx_test) && (has_test)) result[["rfx_preds_test"]] = rfx_preds_test
    class(result) <- "bartmodel"
    
    # Clean up classes with external pointers to C++ data structures
    if (include_mean_forest) rm(forest_model_mean)
    if (include_variance_forest) rm(forest_model_variance)
    rm(forest_dataset_train)
    if (has_test) rm(forest_dataset_test)
    if (has_rfx) rm(rfx_dataset_train, rfx_tracker_train, rfx_model)
    rm(outcome_train)
    rm(rng)
    
    return(result)
}