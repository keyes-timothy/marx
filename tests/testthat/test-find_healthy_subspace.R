healthy_data <-
  matrix(data = runif(n = 300 * 40), nrow = 300, ncol = 40) |>
  dplyr::as_tibble()

# test svd

svd_subspace <-
  healthy_data |>
  marx_find_healthy_subspace(
    method = "svd"
  )

svd_projection <-
  healthy_data |>
  marx_project(healthy_subspace = svd_subspace)

# test nmf

nmf_subspace <-
  healthy_data |>
  marx_find_healthy_subspace(num_metacells = 10L, method = "nmf")

nmf_projection <-
  healthy_data |>
  marx_project(healthy_subspace = nmf_subspace)

# test pma
pma_subspace <-
  healthy_data |>
  marx_find_healthy_subspace(num_metacells = 10L, method = "pma")

pma_projection <-
  healthy_data |>
  marx_project(healthy_subspace = pma_subspace)

# test magnitudes
magnitudes <-
  svd_projection |>
  marx_find_magnitudes()
