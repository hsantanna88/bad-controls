# Extract ATT and SE from ptetools results

Helper to extract the overall ATT estimate and standard error from
either `pte_results` or `pte_emp_boot` objects returned by
[`ptetools::pte_default`](https://rdrr.io/pkg/ptetools/man/pte_default.html).

## Usage

``` r
extract_att(res)
```

## Arguments

- res:

  Result object from
  [`didbc`](https://github.com/hugosantanna/badcontrols/reference/didbc.md)
  or
  [`ptetools::pte_default`](https://rdrr.io/pkg/ptetools/man/pte_default.html)

## Value

A list with components:

- att:

  Overall ATT estimate

- se:

  Standard error

## Examples

``` r
# \donttest{
sim <- simulate_bad_controls(n = 500)
res <- didbc(
  yname = "Y", gname = "G", tname = "period", idname = "id",
  data = sim$data, bad_control_formula = ~X, xformula = ~Z
)
#> Warning: critical value for uniform confidence band is somehow smaller than
#>             critical value for pointwise confidence interval...using pointwise
#>             confidence interal
#> Warning: critical value for uniform confidence band is somehow smaller than
#>             critical value for pointwise confidence interval...using pointwise
#>             confidence interal
#> Warning: critical value for uniform confidence band is somehow smaller than
#>             critical value for pointwise confidence interval...using pointwise
#>             confidence interal
#> Warning: critical value for uniform confidence band is somehow smaller than
#>             critical value for pointwise confidence interval...using pointwise
#>             confidence interal
#> Warning: critical value for uniform confidence band is somehow smaller than
#>             critical value for pointwise confidence interval...using pointwise
#>             confidence interal
extract_att(res)
#> $att
#> [1] 1.57483
#> 
#> $se
#> [1] 0.06555557
#> 
# }
```
