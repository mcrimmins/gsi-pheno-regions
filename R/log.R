# Minimal logger: timestamped lines to the console and to logs/<script>_<profile>_<time>.log

log_init <- function(script, profile) {
  dir.create(here::here("logs"), showWarnings = FALSE)
  f <- here::here("logs", sprintf("%s_%s_%s.log", script, profile,
                                  format(Sys.time(), "%Y%m%d-%H%M%S")))
  options(gsi.log_file = f)
  f
}

log_msg <- function(..., level = "INFO") {
  line <- sprintf("%s [%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level,
                  paste0(...))
  message(line)
  f <- getOption("gsi.log_file")
  if (!is.null(f)) cat(line, "\n", file = f, append = TRUE, sep = "")
  invisible(line)
}

log_warn <- function(...) log_msg(..., level = "WARN")
log_err  <- function(...) log_msg(..., level = "ERROR")
