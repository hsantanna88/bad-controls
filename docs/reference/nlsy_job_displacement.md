# NLSY79 Job Displacement Application Data

A balanced panel of 3,231 NLSY79 respondents observed biennially from
1992 through 2002. The sample applies the paper's positive-earnings
restriction and excludes respondents first treated in 1992, who have no
pre-treatment period within the application window. The outcome is log
earnings, the bad control is an occupation-based wage score, and
treatment is job displacement.

## Usage

``` r
data(nlsy_job_displacement)
```

## Format

A data.frame with 19,386 rows and 8 columns: `id`, `year`,
`log_earnings`, `occ_score`, `group`, `race`, `female`, and
`educ_max_grade`.

## Source

National Longitudinal Survey of Youth 1979 (NLSY79), U.S. Bureau of
Labor Statistics, https://www.nlsinfo.org/investigator/; IPUMS USA,
https://usa.ipums.org/usa/.

## Details

Job displacement is defined as involuntarily leaving a job because of a
layoff/job elimination or plant/company/workplace closure. The
occupation score is time-varying at the individual level because
respondents can change occupations, even though the score is fixed for a
given occupation. Respondents first displaced in 1992 are excluded
because the application requires a pre-treatment period within the
1992–2002 window.

This is a researcher-created subset of public-use NLSY79 data and an
occupation-score merge based on the public-use IPUMS USA 1990 5\\ It is
not an official NLSY79 data release. See the source files in `data-raw/`
for the construction script and provenance.

## Examples

``` r
data(nlsy_job_displacement)
head(nlsy_job_displacement)
#>       id  year log_earnings occ_score group                   race female
#>    <int> <int>        <num>     <num> <int>                 <char> <lgcl>
#> 1:     8  1992     9.798127  2.040221     0 non_black_non_hispanic   TRUE
#> 2:     8  1994     9.928180  2.748872     0 non_black_non_hispanic   TRUE
#> 3:     8  1996    10.085809  2.748872     0 non_black_non_hispanic   TRUE
#> 4:     8  1998     9.998798  2.748872     0 non_black_non_hispanic   TRUE
#> 5:     8  2000    10.218298  2.748872     0 non_black_non_hispanic   TRUE
#> 6:     8  2002    10.357743  2.358675     0 non_black_non_hispanic   TRUE
#>    educ_max_grade
#>             <int>
#> 1:             14
#> 2:             14
#> 3:             14
#> 4:             14
#> 5:             14
#> 6:             14
table(nlsy_job_displacement$group[!duplicated(nlsy_job_displacement$id)])
#> 
#>    0 1994 1996 1998 2000 2002 
#> 2483  209  155  113  103  168 
```
