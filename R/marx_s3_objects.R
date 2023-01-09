# marx_s3_objects.R
# description
#

# constructor

#' TODO?
#'
#' TODO?
#'
#' @param metacells TODO?
#' @param embeddings TODO?
#' @param num_metacells TODO?
#' @param method TODO?
#'
#' @return TODO?
#'
new_marx_subspace <- function(metacells, embeddings, num_metacells, method) {
  healthy_subspace <-
    structure(
      list(
        metacells = metacells,
        embeddings = embeddings,
        num_metacells = num_metacells,
        method = method
      ),
      class = "marx_subspace"
    )

  return(healthy_subspace)
}

# validator --------------------------------------------------------------------

# print method  ----------------------------------------------------------------


#' Print a marx subspace.
#'
#' @param x TODO?
#' @param ... TODO
#'
#' @return TODO?
#'
#' @method print marx_subspace
#'
#' @export
#'
print.marx_subspace <- function(x, ...) {
  print(
    paste0(
      "A marx subspace with ",
      x$num_metacells,
      " metacells constructed using the following method: ",
      toupper(x$method),
      "."
    )
  )

  print("Metacells:")
  print(
    dplyr::as_tibble(x$metacells[,1:x$num_metacells], rownames = "feature")
  )

  print("Embeddings:")
  print(dplyr::as_tibble(x$embeddings[,1:x$num_metacells]))

}


# others?
