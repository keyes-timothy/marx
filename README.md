
<!-- README.md is generated from README.Rmd. Please edit that file -->

# marx

<!-- badges: start -->
<!-- badges: end -->

`marx` is an R package that implements <u>Ma</u>trix factorization and
<u>R</u>esidual E<u>x</u>pression (MARX), an algorithm for comparing
high-dimensional biological measurements to a lower-dimensional
reference linear subspace. Currently, MARX is unpublished and highly
experimental.

## Installation

You can install the development version of marx from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("keyes-timothy/marx")
```

## Example

This is a basic example for using `marx` to decompose an input dataset
into two components: one that aligns with a reference dataset
(represented as a linear subspace) and one that doesn’t align with the
reference.

``` r
library(marx)

# simulate reference data 
# this data is assumed to be an expression matrix in which rows represent 
# observations and columns represent measurements (i.e. genes). 
healthy_data <-
  matrix(data = runif(n = 300 * 40), nrow = 300, ncol = 40) |>
  dplyr::as_tibble()
#> Warning: The `x` argument of `as_tibble.matrix()` must have unique column names if `.name_repair` is omitted as of tibble 2.0.0.
#> Using compatibility `.name_repair`.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.

# simulate new data 
disease_data <- 
  matrix(data = runif(n = 700 * 40, min = 1, max = 2), nrow = 700, ncol = 40) |>
  dplyr::as_tibble()

# build a healthy reference subspace using the SVD
healthy_subspace <- 
  healthy_data |> 
  # other method options include "nmf" and "pma"
  marx_find_healthy_subspace(pca_threshold = 0.9, method = "svd")

# project the disease_data onto the healthy subspace
disease_projection <- 
  disease_data |> 
  marx_project(healthy_subspace = healthy_subspace)

head(disease_projection)
#> # A tibble: 6 × 80
#>   V1_lineage V2_lineage V3_lin…¹ V4_li…² V5_li…³ V6_li…⁴ V7_li…⁵ V8_li…⁶ V9_li…⁷
#>        <dbl>      <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
#> 1       1.71       1.43     1.54    1.62    1.39    1.53    1.52    1.28    1.57
#> 2       1.59       1.58     1.27    1.78    1.43    1.36    1.38    1.64    1.68
#> 3       1.17       1.29     1.78    1.66    1.50    1.85    1.38    1.08    1.51
#> 4       1.33       1.54     1.82    1.29    1.37    1.54    1.50    1.37    1.53
#> 5       1.20       1.57     1.69    1.61    1.52    1.59    1.35    1.58    1.36
#> 6       1.71       1.55     1.32    1.15    1.58    1.61    1.68    1.49    1.59
#> # … with 71 more variables: V10_lineage <dbl>, V11_lineage <dbl>,
#> #   V12_lineage <dbl>, V13_lineage <dbl>, V14_lineage <dbl>, V15_lineage <dbl>,
#> #   V16_lineage <dbl>, V17_lineage <dbl>, V18_lineage <dbl>, V19_lineage <dbl>,
#> #   V20_lineage <dbl>, V21_lineage <dbl>, V22_lineage <dbl>, V23_lineage <dbl>,
#> #   V24_lineage <dbl>, V25_lineage <dbl>, V26_lineage <dbl>, V27_lineage <dbl>,
#> #   V28_lineage <dbl>, V29_lineage <dbl>, V30_lineage <dbl>, V31_lineage <dbl>,
#> #   V32_lineage <dbl>, V33_lineage <dbl>, V34_lineage <dbl>, …
#> # ℹ Use `colnames()` to see all variable names
```

marx can also compute the lineage- and disease-specific component
magnitudes (i.e. the reconstruction and the reconstruction loss for each
observation using the reference subspace):

``` r
disease_projection |> 
  marx_find_magnitudes() |> 
  head()
#> # A tibble: 6 × 3
#>   lineage_magnitude disease_magnitude total_magnitude
#>               <dbl>             <dbl>           <dbl>
#> 1              9.65              1.58            9.78
#> 2              9.69              1.44            9.80
#> 3              9.72              1.74            9.88
#> 4              9.58              1.01            9.64
#> 5              9.53              1.16            9.60
#> 6              9.21              1.44            9.32
```

We can compare the disease-specific component of the healthy data and
the disease data:

``` r
# find healthy projection
healthy_projection <- 
  healthy_data |> 
  marx_project(healthy_subspace = healthy_subspace)

# find healthy magnitudes 
healthy_magnitudes <- 
  healthy_projection |> 
  marx_find_magnitudes()

# find disease magnitudes
disease_magnitudes <- 
  disease_projection |> 
  marx_find_magnitudes()

# make a simple density plot
plot_tibble <- 
  dplyr::bind_rows(
    dplyr::mutate(healthy_magnitudes, sample_type = "healthy"), 
    dplyr::mutate(disease_magnitudes, sample_type = "disease")
  )

plot_tibble |> 
  ggplot2::ggplot(ggplot2::aes(x = disease_magnitude, fill = sample_type)) + 
  ggplot2::geom_density(alpha = 0.5) + 
  ggplot2::theme_bw() + 
  ggplot2::labs(subtitle = "diseased tissue has a larger disease-specific component than healthy tissue")
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />
