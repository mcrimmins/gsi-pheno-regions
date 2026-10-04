# Push notifications via ntfy (https://ntfy.sh or a self-hosted server).
#
# The topic is read from the GSI_NTFY_TOPIC environment variable (set it in ~/.Renviron
# on each machine), never from config.yml: on ntfy.sh anyone who knows the topic can read
# and post to it, so keep it out of git. Optional: GSI_NTFY_SERVER (default from config),
# GSI_NTFY_TOKEN (access token for a protected topic / self-hosted server).
#
# Notifications never stop the pipeline: failures to send are reported and ignored.

`%||%` <- function(a, b) if (is.null(a)) b else a

notify_init <- function(cfg) {
  n <- cfg$notify %||% list()
  topic <- Sys.getenv("GSI_NTFY_TOPIC", "")
  server <- Sys.getenv("GSI_NTFY_SERVER", n$server %||% "https://ntfy.sh")
  st <- list(
    on = isTRUE(n$enabled) && nzchar(topic),
    url = paste0(sub("/+$", "", server), "/", topic),
    token = Sys.getenv("GSI_NTFY_TOKEN", ""),
    tag = sprintf("%s/%s", Sys.info()[["nodename"]], cfg$profile),
    downloads = isTRUE(n$download_updates),
    max_errors = as.integer(n$max_errors %||% 10L),
    n_err = 0L
  )
  options(gsi.notify = st)
  if (isTRUE(n$enabled) && !st$on) {
    log_msg("ntfy: GSI_NTFY_TOPIC not set; notifications off")
  } else if (st$on) {
    log_msg("ntfy: notifications on via ", sub("/+$", "", server))
  }
  invisible(st)
}

.ascii <- function(x) iconv(x, to = "ASCII//TRANSLIT", sub = "?")

# priority: 1 min, 2 low, 3 default, 4 high, 5 urgent. tags: ntfy emoji short codes.
notify <- function(title, msg, priority = 3, tags = NULL) {
  st <- getOption("gsi.notify")
  if (is.null(st) || !isTRUE(st$on)) return(invisible(FALSE))
  req <- httr2::request(st$url) |>
    httr2::req_headers(Title = .ascii(sprintf("[%s] %s", st$tag, title)),
                       Priority = as.character(priority)) |>
    httr2::req_body_raw(enc2utf8(msg), type = "text/plain; charset=utf-8") |>
    httr2::req_timeout(15) |>
    httr2::req_error(is_error = function(resp) FALSE)
  if (length(tags)) req <- httr2::req_headers(req, Tags = paste(tags, collapse = ","))
  if (nzchar(st$token)) req <- httr2::req_auth_bearer_token(req, st$token)
  ok <- tryCatch(httr2::resp_status(httr2::req_perform(req)) < 300,
                 error = function(e) FALSE)
  if (!ok) message("ntfy: notification failed (continuing)")
  invisible(ok)
}

# Errors are capped per run so an outage can't flood the phone.
notify_error <- function(title, msg, priority = 4) {
  st <- getOption("gsi.notify")
  if (is.null(st) || !isTRUE(st$on)) return(invisible(FALSE))
  st$n_err <- st$n_err + 1L
  options(gsi.notify = st)
  if (st$n_err <= st$max_errors) {
    notify(title, msg, priority = priority, tags = "warning")
  } else if (st$n_err == st$max_errors + 1L) {
    notify("further errors muted",
           sprintf("More than %d errors this run; see the log.", st$max_errors),
           priority = 4, tags = "mute")
  }
}

notify_downloads_on <- function() isTRUE(getOption("gsi.notify")$downloads)

fmt_dur <- function(secs) {
  if (secs < 90) sprintf("%.0f s", secs) else if (secs < 5400) sprintf("%.0f min", secs / 60)
  else sprintf("%.1f h", secs / 3600)
}
