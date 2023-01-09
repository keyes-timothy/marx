# find_healthy_subspace.R
#
#

# Find reference subspace ------------------------------------------------------

#' Find a low-dimensional subspace representing an approximate basis for
#' healthy cytometry data using the SVD
#'
#' @param healthy_data A data.frame or matrix of dimension (n x p) in which rows
#' represent cells and columns represent single-cell measurements
#' (i.e. CyTOF antigens).
#' Only cells taken from healthy (or reference) samples should be included.
#'
#' @param pca_threshold a value between 0 and 1 representing the percentage of
#' variance in the original data that should be preserved during PCA transformation.
#' This will determine the number k of PCs (i.e. metacells) in the result.
#'
#' @param return_recipe A boolean value indicating if the recipe for the pca
#' should be returned (TRUE) (or if the projection itself should be returned; FALSE; the
#' default). Not currently used.
#'
#' @return A `marx_subspace` S3 object with 4 slots.
#' `metacells` is a (p x k) matrix in which the
#' columns represent the k metacells and the p rows represent the p features in the
#' input data. Each entry represents the weighting of a feature on each metacell.
#' `embeddings` is a (n x k) matrix giving the projection of each of the n
#' input cells onto the k metacells (these values are often called the principal
#' components, i.e. the first column contains the scores for each of the n input
#' cells on the first principal component).
#' `num_metacells` is an integer of length 1 indicating how many PCs/metacells
#' are needed to preserve the proportion of variance specified by `pca_threshold`.
#' `method` is a character vector of length 1 indicating which matrix factorization
#' method was used (in this case, "svd").
#'
#'
#'
marx_find_healthy_subspace_svd <-
  function(
    healthy_data,
    pca_threshold = 0.95,
    return_recipe = FALSE
  ) {

    # compute the svd of the input data matrix
    healthy_svd <- svd(x = as.matrix(healthy_data))

    # extract the right singular vectors of the svd, which can be interpreted
    # as a num_features x num_metacells weighting of each input feature on each
    # output metacell
    eigencells <- healthy_svd$v
    colnames(eigencells) <- paste0("metacell_", 1:ncol(eigencells))
    rownames(eigencells) <- colnames(healthy_data)

    # extract the pca embeddings for each cell in the input dataset in the
    # new subspace
    pca_embeddings <- healthy_svd$u %*% diag(healthy_svd$d)
    colnames(pca_embeddings) <- paste0("PC", 1:ncol(eigencells))

    # estimate standard deviations using nrow(healthy_data) - 1 as the
    # denominator, as in prcomp()
    pca_standard_deviations <-
      healthy_svd$d / sqrt(max(1, nrow(healthy_data) - 1))
    pca_variances <- pca_standard_deviations^2

    # calculate how many eigencells are needed to preserve the specified
    # amount of variance in the input data
    cumulative_prop_variance <-
      (pca_variances / sum(pca_variances)) |>
      cumsum()
    num_eigencells <- sum(cumulative_prop_variance < pca_threshold) + 1

    # return result
    healthy_subspace <-
      new_marx_subspace(
        metacells = eigencells,
        embeddings = pca_embeddings,
        num_metacells = num_eigencells,
        method = "svd"
      )

    return(healthy_subspace)
  }

#' Find a low-dimensional subspace representing an approximate basis for
#' healthy cytometry data using non-negative matrix factorization (NMF).
#'
#' @param healthy_data A data.frame or matrix of dimension (n x p) in which rows
#' represent cells and columns represent single-cell measurements
#' (i.e. CyTOF antigens).
#' Only cells taken from healthy (or reference) samples should be included.
#'
#' @param num_metacells The number of metacells (i.e. the number of NMF
#' components) to use in constructing the healthy subspace.
#'
#' @return A `marx_subspace` S3 object with 4 slots.
#' `metacells` is a (p x k) matrix in which the
#' columns represent the k metacells and the p rows represent the p features in the
#' input data. Each entry represents the weighting of a feature on each metacell.
#' `embeddings` is a (n x k) matrix giving the projection of each of the n
#' input cells onto the k metacells (these values are often called NMF embeddings,
#' i.e. the first column contains the scores for each of the n input
#' cells on the first NMF factor).
#' `num_metacells` is an integer of length 1 indicating how many metacells
#' were specified in the original function call.
#' `method` is a character vector of length 1 indicating which matrix factorization
#' method was used (in this case, "nmf").
#'
#' @importFrom RcppML nmf
#'
#' @export
#'
marx_find_healthy_subspace_nmf <-
  function(
    healthy_data,
    num_metacells = 10L
  ) {
    # nmf_recipe <-
    #   recipes::recipe(healthy_data) |>
    #   recipes::step_nnmf_sparse(recipes::all_numeric(), num_comp = num_metacells) |>
    #   recipes::prep()

    nmf_result <-
      RcppML::nmf(
        A = t(as.matrix(healthy_data)),
        k = num_metacells,
        verbose = FALSE
      )

    # extact the nmf embeddings for each input cell, i.e. the projection
    # of each cell onto each of the num_metacells metacells
    # factor_matrix <-
    #   nmf_recipe |>
    #   recipes::juice() |>
    #   as.matrix()
    nmf_embeddings <-
      nmf_result$h |>
      t()
    colnames(nmf_embeddings) <- paste0("metacell_", 1:ncol(nmf_embeddings))

    # extract the metacells
    # in this case, the metacells are the same as the feature factor
    # matrix from RCppML::nmf (aka w or h, depending on whether or not you
    # transpose the input data), and they form a basis
    # for the healthy subspace. Each col can be interpreted as a cell state.
    #metacells <- nmf_recipe$steps[[1]]$res$w
    metacells <- nmf_result$w
    rownames(metacells) <- colnames(healthy_data)
    colnames(metacells) <- paste0("metacell_", 1:ncol(metacells))

    # return result
    healthy_subspace <-
      new_marx_subspace(
        metacells = metacells,
        embeddings = nmf_embeddings,
        num_metacells = num_metacells,
        method = "nmf"
      )

    return(healthy_subspace)
  }



#' Find a low-dimensional subspace representing an approximate basis for
#' healthy cytometry data using penalized matrix analysis (PMA).
#'
#' @param healthy_data A data.frame or matrix of dimension (n x p) in which rows
#' represent cells and columns represent single-cell measurements
#' (i.e. CyTOF antigens).
#' Only cells taken from healthy (or reference) samples should be included.
#'
#' @param num_metacells The number of metacells (i.e. the number of PMA
#' components) to use in constructing the healthy subspace.
#'
#' @param sumabsv A numeric vector of length 1 indicating how sparse the
#' decomposition should be (it is equal to the sum of the absolute value of the
#' elements of v, PMA's sparse principal components).
#'
#' @return A `marx_subspace` S3 object with 4 slots.
#' `metacells` is a (p x k) matrix in which the
#' columns represent the k metacells and the p rows represent the p features in the
#' input data. Each entry represents the weighting of a feature on each metacell.
#' `embeddings` is a (n x k) matrix giving the projection of each of the n
#' input cells onto the k metacells (these values are often called sparse
#' principal component embeddings,
#' i.e. the first column contains the scores for each of the n input
#' cells on the first sparse principal component).
#' `num_metacells` is an integer of length 1 indicating how many metacells
#' were specified in the original function call.
#' `method` is a character vector of length 1 indicating which matrix factorization
#' method was used (in this case, "pma").
#'
#' @export
#'
#' @importFrom PMA SPC
#'
marx_find_healthy_subspace_pma <-
  function(
    healthy_data = NULL,
    num_metacells = 10,
    sumabsv
  ) {
    # keep track of the antigen names in the input tibble
    col_names <-
      healthy_data |>
      colnames()

    v_start <- svd(x = as.matrix(healthy_data))$v

    if (!missing(sumabsv)) {
      best_sumabsv <- sumabsv
    } else {
      best_sumabsv <-
        # PMA::SPC.cv(
        #   x = as.matrix(healthy_data),
        #   sumabsvs = seq(1, sqrt(ncol(healthy_data)), length.out = 5),
        #   v = v_start,
        #   vpos = TRUE
        # )$bestsumabsv1se
        sqrt(ncol(healthy_data)) / 2
    }

    pma_result <-
      PMA::SPC(
        x = as.matrix(healthy_data),
        sumabsv = best_sumabsv,
        K = num_metacells,
        cnames = colnames(healthy_data),
        v = v_start,
        vpos = TRUE,
        compute.pve = FALSE
      )

    pma_embeddings <- pma_result$u

    metacells <- marx_extract_metacells_pma(pma_result)

    # return result
    healthy_subspace <-
      new_marx_subspace(
        metacells = metacells,
        embeddings = pma_embeddings,
        num_metacells = num_metacells,
        method = "pma"
      )

    return(healthy_subspace)
  }


marx_extract_metacells_pma <-
  function(marx_pma_result) {
    metacells <- marx_pma_result$v
    colnames(metacells) <- paste0("PMA", 1:ncol(metacells))
    rownames(metacells) <- marx_pma_result$cnames
    return(metacells)
  }


#' Find a low-dimensional subspace representing an approximate basis for
#' healthy single-cell data using matrix factorization
#'
#' @param healthy_data A data.frame or matrix of dimension (n x p) in which rows
#' represent cells and columns represent single-cell measurements
#' (i.e. CyTOF antigens).
#' Only cells taken from healthy (or reference) samples should be included.
#' @param ... Additional arguments to pass to the specific method.
#' @param method A character vector of length 1 indicating which matrix
#' factorization method should be used. Valid options are "svd", "nmf", and "pma".
#'
#' @return A `marx_subspace` S3 object with 4 slots.
#' `metacells` is a (p x k) matrix in which the
#' columns represent the k metacells and the p rows represent the p features in the
#' input data. Each entry represents the weighting of a feature on each metacell.
#' `embeddings` is a (n x k) matrix giving the projection of each of the n
#' input cells onto the k metacells.
#' `num_metacells` is an integer of length 1 indicating how many metacells
#' were requested for the subspace representation.
#' `method` is a character vector of length 1 indicating which matrix factorization
#' method was used.
#'
#' @export
#'
#'
marx_find_healthy_subspace <-
  function(
    healthy_data,
    ...,
    method = c("svd", "nmf", "pma")
  ) {

    if (method == "svd") {
      result <- marx_find_healthy_subspace_svd(healthy_data = healthy_data, ...)
    }

    else if (method == "nmf") {
      result <- marx_find_healthy_subspace_nmf(healthy_data = healthy_data, ...)
    }

    else if (method == "pma") {
      result = marx_find_healthy_subspace_pma(healthy_data = healthy_data, ...)
    }

    else {
      stop("Not a valid method.")
    }

    return(result)
  }


# project onto reference subspace ----------------------------------------------


#' Project single-cell measurements onto a PCA-defined healthy subspace
#'
#' @param new_data A tibble in which rows represent cells and columns represent
#' antigens measured for each cell. Only the columns in `new_data` corresponding
#' to antigens used to build the healthy subspace will be used in the projection.
#' All others will be ignored.
#'
#' @param healthy_subspace A matrix output of the function `marx_find_healthy_subspace`
#'
#' @return A tibble with an equal number of rows as `new_data` and twice the number
#' of columns. Each row will represent a cell from `new_data`, and each column
#' will represent either the healthy or disease-specific component of the
#' healthy subspace projection.
#'
#' NEEDS TO CHECK THAT THE HEALTHY BASIS IS COMPATIBLE WITH THE new_data IF
#' CLUSTER_COL IS SPECIFIED
#'
#' @importFrom dplyr all_of
#' @importFrom dplyr as_tibble
#' @importFrom dplyr mutate
#' @importFrom dplyr rename_with
#' @importFrom dplyr select
#'
#' @importFrom stats lm
#'
#' @export
#'
marx_project <-
  function(
    new_data,
    healthy_subspace
  ) {

    # appease R CMD check
    cell_id <- NULL

    if (inherits(healthy_subspace, what = "matrix")) {
      antigens <- row.names(healthy_subspace)
    } else {
      antigens <- row.names(healthy_subspace$metacells)
    }

    # check that new_data and healthy_subspace are compatible
    if (!all(antigens %in% colnames(new_data))) {
      stop("the columns of new_data do not contain all the antigens used
           to create the healthy subspace")
    }

    # reformat the input data to be projected as a [p x n] matrix
    cell_matrix <-
      new_data |>
      # give each cell an id so that the output tibble has rows in the
      # same order as the input tibble
      dplyr::mutate(cell_id = 1:nrow(new_data)) |>
      dplyr::select(cell_id, dplyr::all_of(antigens)) |>
      # transpose matrix and name rows/columns
      cell_matrix_transformation()

    # only project onto the requested number of metacells
    loadings <- healthy_subspace$metacells[, 1:healthy_subspace$num_metacells]

    # to perform the projection,
    # fit one linear model for each cell (column) in the cell matrix
    # without an intercept
    lm_result <- stats::lm(cell_matrix ~ loadings + 0)

    # find healthy component
    healthy_component <-
      lm_result$fitted.values |>
      t() |>
      dplyr::as_tibble() |>
      dplyr::rename_with(.fn = function(x) paste0(x, "_lineage"))

    # find disease component
    disease_component <-
      lm_result$residuals |>
      t() |>
      dplyr::as_tibble() |>
      dplyr::rename_with(.fn = function(x) paste0(x, "_disease"))

    # combine components into the final result

    result <- dplyr::bind_cols(healthy_component, disease_component)

    # return result
    return(result)
  }


# will return the cell_matrix such that the column names correspond to each
# cell's cell_id in the original input tibble.
cell_matrix_transformation <- function(cell_matrix) {

  result <-
    cell_matrix |>
    t()

  colnames(result) <- result[1,]
  result <- result[2:nrow(result), ]

  return(result)
}


