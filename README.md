# Book for scrapper

This repository contains a Bioconductor package to deploy the [**scrapper**](https://github.com/libscran/scrapper) book.
Install all relevant packages used in the book with:

```r
BiocManager::install(remotes::local_package_deps(dependencies=TRUE))
```

Building the book can be done by either running the usual **bookdown** invocation in `inst/book`:

```r
bookdown::render_book("index.Rmd")
```

Or by `R CMD build` on the package itself to compile the book during the vignette build process.
