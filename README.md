# KMerEnrich 2.0 — development branch

> **Work in progress — not a release.** This branch contains the in-development
> version 2.0 of KmerEnrich. It is not tagged, may change at any time, and should
> not be cited.
>
> The released and citable version is **[v0.1.0](https://github.com/Miloumille/KmerEnrich/releases/tag/v0.1.0)**
> on the `main` branch.

## What changes in 2.0

Version 2.0 is a rewrite focused on designing semi-specific SMART-9N primers that
enrich viral RNA over a host/vector background, using k-mer statistics,
thermodynamic filtering and integer-linear-programming coverage optimisation, with
an interactive HTML report.

## Install (development version)

```r
devtools::install_github("Miloumille/KmerEnrich@dev")
```

See `vignettes/getting-started.Rmd` and the example scripts in `scripts/`.

## License

MIT. See [LICENSE](LICENSE).
