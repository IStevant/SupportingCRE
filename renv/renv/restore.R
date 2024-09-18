source(".Rprofile")
# Install all the missing R packages from the renv.lock file
renv::restore()