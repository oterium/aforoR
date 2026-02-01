#' Select Discriminative Points for Functional Data Classification
#'
#' Implements the Hall & Bathia (2012) algorithm for selecting a subset of points
#' from functional data (e.g., wavelet coefficients or contours) that optimize
#' classification accuracy.
#'
#' @param data An object of class \code{fdata} or a matrix where each row is an observation.
#' @param grouping A factor specifying the class labels for each observation.
#' @param method Classification method to use. One of "lda", "qda", "NaiveBayes",
#'               "logistic", "knn", "svm".
#' @param cv Logical. If TRUE, use cross-validation to estimate error.
#' @param p Numeric. Stopping criterion. Stop when improvement is less than p * error.
#' @param delta_t Numeric. Minimum distance between selected points (as a fraction of total points).
#' @param mean_error Metric for error. Either "mtotal" (overall error) or "mgroups" (mean error per group).
#' @param parallel Logical. If TRUE, use parallel processing (supported on Windows via \code{parLapply} or Unix via \code{mclapply}).
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{selected_indices}: Indices of the selected points.
#'     \item \code{error_path}: Vector of error rates at each iteration.
#'     \item \code{confusion_matrix}: Final confusion matrix using cross-validation.
#'   }
#' @export
#' @importFrom stats predict na.omit
#' @examples
#' \dontrun{
#' # Example usage with dummy data
#' library(fda.usc)
#' data(phoneme)
#' res <- select_points_hall(phoneme$data[1:100, ], phoneme$class[1:100], method = "lda")
#' }
select_points_hall <- function(data, grouping, method = "lda", cv = TRUE,
                               p = 0.01, delta_t = 0.01, mean_error = "mgroups",
                               parallel = FALSE) {
    # Ensure fda.usc is available if data is fdata
    if (inherits(data, "fdata")) {
        if (!requireNamespace("fda.usc", quietly = TRUE)) {
            stop("Package 'fda.usc' is required for fdata objects. Please install it.")
        }
        x_data <- data$data
    } else {
        x_data <- as.matrix(data)
    }

    grouping <- as.factor(grouping)
    n_points <- ncol(x_data)
    delta_idx <- max(1, floor(delta_t * n_points))

    # helper for error calculation
    calc_ferror <- function(groups, pred, metric) {
        tab <- table(groups, pred)
        if (metric == "mtotal") {
            return(1 - sum(diag(tab)) / sum(tab))
        } else {
            # mgroups: mean error of groups
            diag_vals <- diag(tab)
            row_sums <- rowSums(tab)
            # avoid division by zero
            group_errors <- 1 - (diag_vals / ifelse(row_sums == 0, 1, row_sums))
            return(mean(group_errors))
        }
    }

    # Wrapper for classification methods
    run_discriminant <- function(x_train, groups, method, cv_opt) {
        # Ensure necessary packages are available
        pkg <- switch(method,
            "lda" = "MASS",
            "qda" = "MASS",
            "NaiveBayes" = "klaR",
            "logistic" = "nnet",
            "knn" = "class",
            "svm" = "e1071"
        )

        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop(sprintf("Package '%s' is required for method '%s'.", pkg, method))
        }

        if (method == "lda") {
            if (cv_opt) {
                fit <- MASS::lda(x_train, grouping = groups, CV = TRUE)
                return(fit$class)
            } else {
                fit <- MASS::lda(x_train, grouping = groups)
                return(predict(fit, x_train)$class)
            }
        } else if (method == "qda") {
            if (cv_opt) {
                fit <- MASS::qda(x_train, grouping = groups, CV = TRUE)
                return(fit$class)
            } else {
                fit <- MASS::qda(x_train, grouping = groups)
                return(predict(fit, x_train)$class)
            }
        } else if (method == "NaiveBayes") {
            # Naive Bayes CV manual implementation (as in original)
            if (cv_opt) {
                preds <- sapply(1:length(groups), function(i) {
                    fit <- klaR::NaiveBayes(x_train[-i, , drop = FALSE], grouping = groups[-i])
                    as.character(predict(fit, x_train[i, , drop = FALSE])$class)
                })
                return(factor(preds, levels = levels(groups)))
            } else {
                fit <- klaR::NaiveBayes(x_train, grouping = groups)
                return(predict(fit, x_train)$class)
            }
        } else if (method == "logistic") {
            if (cv_opt) {
                preds <- sapply(1:length(groups), function(i) {
                    fit <- nnet::multinom(groups[-i] ~ ., data = as.data.frame(x_train[-i, , drop = FALSE]), trace = FALSE)
                    as.character(predict(fit, as.data.frame(x_train[i, , drop = FALSE])))
                })
                return(factor(preds, levels = levels(groups)))
            } else {
                fit <- nnet::multinom(groups ~ ., data = as.data.frame(x_train), trace = FALSE)
                return(predict(fit, as.data.frame(x_train)))
            }
        } else if (method == "knn") {
            k_val <- length(levels(groups))
            if (cv_opt) {
                return(class::knn.cv(train = x_train, cl = groups, k = k_val))
            } else {
                return(class::knn(train = x_train, test = x_train, cl = groups, k = k_val))
            }
        } else if (method == "svm") {
            # Simple SVM without tuning for performance during selection
            if (cv_opt) {
                preds <- sapply(1:length(groups), function(i) {
                    fit <- e1071::svm(x_train[-i, , drop = FALSE], groups[-i])
                    as.character(predict(fit, x_train[i, , drop = FALSE]))
                })
                return(factor(preds, levels = levels(groups)))
            } else {
                fit <- e1071::svm(x_train, groups)
                return(predict(fit, x_train))
            }
        }
    }

    find_best_next_point <- function(available_indices, current_indices) {
        eval_point <- function(idx) {
            test_indices <- c(current_indices, idx)
            preds <- run_discriminant(x_data[, test_indices, drop = FALSE], grouping, method, cv)
            calc_ferror(grouping, preds, mean_error)
        }

        if (parallel) {
            if (.Platform$OS.type == "unix") {
                errors <- parallel::mclapply(available_indices, eval_point)
            } else {
                # Windows parallel
                cl <- parallel::makeCluster(parallel::detectCores() - 1)
                on.exit(parallel::stopCluster(cl))
                parallel::clusterExport(cl, varlist = c("x_data", "grouping", "method", "cv", "mean_error", "run_discriminant", "calc_ferror"), envir = environment())
                errors <- parallel::parLapply(cl, available_indices, eval_point)
            }
        } else {
            errors <- lapply(available_indices, eval_point)
        }

        errors <- unlist(errors)
        best_idx <- available_indices[which.min(errors)]
        return(list(error = min(errors), index = best_idx))
    }

    # Sequential phase
    M <- 1:n_points
    selected_tn <- c()
    error_history <- c()

    # r=1
    step1 <- find_best_next_point(M, c())
    error_history <- c(error_history, step1$error)
    selected_tn <- c(selected_tn, step1$index)

    # Iterative selection
    while (TRUE) {
        # Filter available indices based on delta_t constraint
        last_point <- selected_tn[length(selected_tn)]
        M_available <- M[!(abs(M - last_point) < (2 * delta_idx))]
        # Also exclude already selected points if they weren't excluded by delta_t
        M_available <- setdiff(M_available, selected_tn)

        if (length(M_available) == 0) break

        step <- find_best_next_point(M_available, selected_tn)

        # Stopping criterion
        improvement <- (error_history[length(error_history)] - step$error)
        if (improvement < p * error_history[length(error_history)]) {
            break
        }

        error_history <- c(error_history, step$error)
        selected_tn <- c(selected_tn, step$index)

        # QDA constraint: can't have more points than samples in smallest group - 2
        if (method == "qda") {
            min_group_size <- min(table(grouping))
            if (length(selected_tn) >= min_group_size - 2) break
        }

        # Limit maximum points to avoid overfitting and computational cost
        if (length(selected_tn) >= 10) break
    }

    # Refinement Phase (simplified and modernized)
    if (length(selected_tn) >= 2 && length(selected_tn) <= 4) {
        bufer <- if (length(selected_tn) <= 3) 10 else 5

        # Create search ranges around selected points
        ranges <- lapply(selected_tn, function(idx) {
            low <- max(1, idx - 2 * delta_idx * bufer)
            high <- min(n_points, idx + 2 * delta_idx * bufer)
            seq(low, high, by = max(1, 2 * delta_idx))
        })

        # Generate combinations using expand.grid
        grid <- do.call(expand.grid, ranges)

        # Function to evaluate a combination
        eval_comb <- function(row_idx) {
            indices <- as.numeric(grid[row_idx, ])
            preds <- run_discriminant(x_data[, indices, drop = FALSE], grouping, method, cv)
            calc_ferror(grouping, preds, mean_error)
        }

        if (parallel) {
            # Re-use or create cluster for grid search if necessary
            # For simplicity, we use same logic as above
            if (.Platform$OS.type == "unix") {
                grid_errors <- parallel::mclapply(1:nrow(grid), eval_comb)
            } else {
                cl <- parallel::makeCluster(parallel::detectCores() - 1)
                on.exit(parallel::stopCluster(cl))
                parallel::clusterExport(cl, varlist = c("x_data", "grouping", "method", "cv", "mean_error", "run_discriminant", "calc_ferror", "grid"), envir = environment())
                grid_errors <- parallel::parLapply(cl, 1:nrow(grid), eval_comb)
            }
        } else {
            grid_errors <- lapply(1:nrow(grid), eval_comb)
        }

        grid_errors <- unlist(grid_errors)
        best_row <- which.min(grid_errors)

        if (grid_errors[best_row] < error_history[length(error_history)]) {
            error_history[length(error_history)] <- grid_errors[best_row]
            selected_tn <- as.numeric(grid[best_row, ])
        }
    }

    # Final cross-validation confusion matrix
    final_preds <- run_discriminant(x_data[, selected_tn, drop = FALSE], grouping, method, cv_opt = TRUE)
    conf_matrix <- table(groups = grouping, predicted = final_preds)

    return(list(
        selected_indices = selected_tn,
        error_path = error_history,
        confusion_matrix = conf_matrix
    ))
}
