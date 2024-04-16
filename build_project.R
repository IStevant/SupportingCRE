if(!require(fs)) {
	install.packages("fs")
}

if (!require("pacman")) {
	install.packages("pacman")
}

suppressPackageStartupMessages(library(fs))

folders <- c(
	"scripts",
	"data",
	"processed_data",
	"results",
	"graphs"
)

fs::dir_create(
	path = folders
)

if(!require(renv)) {
	install.packages("renv")
}

renv::init()
renv::snapshot()