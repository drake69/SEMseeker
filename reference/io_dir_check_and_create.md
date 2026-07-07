# Create a directory path, building any missing intermediate directories

Splits `baseFolder` into its path components, appends `subFolders`, and
creates each directory level that does not yet exist. Equivalent to
`dir.create(path, recursive = TRUE)` but also returns the final
normalised absolute path.

## Usage

``` r
io_dir_check_and_create(baseFolder, subFolders)
```

## Arguments

- baseFolder:

  Character scalar: root directory path (need not exist yet).

- subFolders:

  Character vector: one or more subdirectory names to append below
  `baseFolder`. Each element becomes one level of the hierarchy.

## Value

Character scalar: the normalised absolute path of the deepest directory
created (or already existing).
