# Build a polynomial regression formula with covariate interactions

Pure helper: no I/O, no side effects. Constructs a formula of the form
`y ~ I(x^1) + I(x^2) + ... + I(x^degree) + I(x^1):cov1 + I(x^1):cov2 + ...`

## Usage

``` r
assoc_polynomial_formula_build(
  dependent_variable,
  independent_variable,
  degree,
  covariates
)
```

## Arguments

- dependent_variable:

  Name of the response column (string).

- independent_variable:

  Name of the predictor column (string).

- degree:

  Polynomial degree (positive integer).

- covariates:

  Character vector of covariate column names (may be empty).

## Value

A `formula` object.
