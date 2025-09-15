all: compiled

compiled: 
	"${R_HOME}/bin/R" -e "setwd('book'); bookdown::render_book('index.Rmd'); dir.create('../../inst/doc', recursive=TRUE); file.rename('docs', '../../inst/doc/book')"
