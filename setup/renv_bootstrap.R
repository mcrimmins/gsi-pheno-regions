# One-time renv setup on the development machine (run from the project root):
#   Rscript setup/renv_bootstrap.R
# Other machines (the P720): open R in the project root and run renv::restore().
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", repos = "https://cloud.r-project.org")
if (!file.exists("renv.lock")) {
  renv::init(bare = TRUE, restart = FALSE)
}
renv::install(c("config", "here", "httr2", "terra", "future", "furrr"))
renv::snapshot(type = "implicit", prompt = FALSE)
cat("renv ready:", renv::paths$library(), "\n")
