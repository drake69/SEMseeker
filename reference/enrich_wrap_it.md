# Wrap long strings to a fixed width

Wraps each string in `x` to at most `len` characters per line,
collapsing the result with newlines. Useful for plot axis labels.

## Usage

``` r
enrich_wrap_it(x, len)
```

## Arguments

- x:

  Character vector of strings to wrap.

- len:

  Integer maximum line width in characters.

## Value

Character vector of the same length as `x`, with long strings broken
into multiple lines separated by `\n`.
