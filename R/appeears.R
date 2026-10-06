# Minimal AppEEARS API client (https://appeears.earthdatacloud.nasa.gov/api).
#
# Credentials: NASA Earthdata login read from a netrc file, never from code or config:
#   machine urs.earthdata.nasa.gov login <user> password <password>
# Looked up in $GSI_NETRC, then ~/.netrc, ~/_netrc, %USERPROFILE%/.netrc, %USERPROFILE%/_netrc.
# Keep it private (chmod 600 on Linux). The token from /login is valid for 48 h.

ae_netrc_creds <- function(host = "urs.earthdata.nasa.gov") {
  cands <- c(Sys.getenv("GSI_NETRC"), path.expand(c("~/.netrc", "~/_netrc")),
             file.path(Sys.getenv("USERPROFILE"), c(".netrc", "_netrc")))
  cands <- cands[nzchar(cands) & file.exists(cands)]
  for (f in cands) {
    tok <- scan(f, what = "", quiet = TRUE)
    i <- which(tok == "machine" & c(tok[-1], "") == host)
    if (length(i)) {
      rest <- tok[(i[1] + 2):length(tok)]
      stop_at <- which(rest == "machine")
      if (length(stop_at)) rest <- rest[seq_len(stop_at[1] - 1)]
      login <- rest[which(rest == "login") + 1]
      pass <- rest[which(rest == "password") + 1]
      if (length(login) && length(pass)) return(list(user = login, password = pass, file = f))
    }
  }
  stop("No Earthdata credentials for ", host, " in a netrc file. Create ~/.netrc (Linux) ",
       "or %USERPROFILE%/_netrc (Windows) with: machine ", host,
       " login <user> password <password>, or set GSI_NETRC to its path.")
}

ae_req <- function(cfg, path, token = NULL) {
  r <- httr2::request(paste0(sub("/+$", "", cfg$cutb$api), "/", path)) |>
    httr2::req_user_agent(cfg$prism$user_agent) |>
    httr2::req_timeout(300) |>
    httr2::req_options(http_version = 2L)
  if (!is.null(token)) r <- httr2::req_auth_bearer_token(r, token)
  r
}

ae_login <- function(cfg) {
  cr <- ae_netrc_creds()
  resp <- ae_req(cfg, "login") |>
    httr2::req_auth_basic(cr$user, cr$password) |>
    httr2::req_method("POST") |>
    httr2::req_perform()
  httr2::resp_body_json(resp)$token
}

ae_logout <- function(cfg, token) {
  try(ae_req(cfg, "logout", token) |> httr2::req_method("POST") |>
        httr2::req_error(is_error = function(r) FALSE) |> httr2::req_perform(), silent = TRUE)
  invisible(NULL)
}

# Layer metadata for a product (public endpoint). Named list: layer -> list(...).
ae_product_layers <- function(cfg, product) {
  resp <- ae_req(cfg, paste0("product/", product)) |>
    httr2::req_error(is_error = function(r) FALSE) |>
    httr2::req_perform()
  if (httr2::resp_status(resp) == 404) {
    stop("Product '", product, "' is not available in AppEEARS. List what is with: ",
         "vapply(httr2::resp_body_json(httr2::req_perform(httr2::request('",
         sub("/+$", "", cfg$cutb$api), "/product'))), function(x) x$ProductAndVersion, '')")
  }
  if (httr2::resp_status(resp) >= 400) stop("AppEEARS product query failed: HTTP ", httr2::resp_status(resp))
  httr2::resp_body_json(resp)
}

# Resolve the configured short names (regex patterns) to exact layer names; stop clearly
# if a pattern matches zero or several layers.
ae_resolve_layers <- function(cfg, product) {
  meta <- ae_product_layers(cfg, product)
  pats <- unlist(cfg$cutb$layers)
  out <- vapply(names(pats), function(nm) {
    hit <- grep(pats[[nm]], names(meta), value = TRUE)
    if (length(hit) != 1) {
      stop(sprintf("%s: pattern '%s' (%s) matched %d layers. Available: %s", product,
                   pats[[nm]], nm, length(hit), paste(names(meta), collapse = ", ")))
    }
    hit
  }, "")
  list(layers = out, meta = meta[out])
}

ae_bbox_geojson <- function(b) {
  ring <- list(c(b[["xmin"]], b[["ymin"]]), c(b[["xmax"]], b[["ymin"]]),
               c(b[["xmax"]], b[["ymax"]]), c(b[["xmin"]], b[["ymax"]]),
               c(b[["xmin"]], b[["ymin"]]))
  list(type = "FeatureCollection",
       features = list(list(type = "Feature", properties = setNames(list(), character()),
                            geometry = list(type = "Polygon", coordinates = list(ring)))))
}

# Submit one area task: one product, all its layers, one calendar year. Returns task_id.
ae_submit_year <- function(cfg, token, product, layers, year, bbox, task_name) {
  body <- list(
    task_type = "area",
    task_name = task_name,
    params = list(
      dates = list(list(startDate = sprintf("01-01-%d", year),
                        endDate = sprintf("12-31-%d", year))),
      layers = lapply(unname(layers), function(l) list(product = product, layer = l)),
      output = list(format = list(type = "geotiff"), projection = "geographic"),
      geo = ae_bbox_geojson(bbox)
    )
  )
  resp <- ae_req(cfg, "task", token) |>
    httr2::req_body_json(body, auto_unbox = TRUE) |>
    httr2::req_perform()
  httr2::resp_body_json(resp)$task_id
}

ae_task_status <- function(cfg, token, task_id) {
  httr2::resp_body_json(httr2::req_perform(ae_req(cfg, paste0("task/", task_id), token)))$status
}

ae_bundle_files <- function(cfg, token, task_id) {
  b <- httr2::resp_body_json(httr2::req_perform(ae_req(cfg, paste0("bundle/", task_id), token)))
  if (!length(b$files)) {
    return(data.frame(file_id = character(), file_name = character(), file_size = numeric(),
                      sha256 = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, lapply(b$files, function(f) data.frame(
    file_id = f$file_id, file_name = basename(f$file_name),
    file_size = as.numeric(f$file_size %||% NA),
    sha256 = f$sha256 %||% NA_character_, stringsAsFactors = FALSE)))
}

# Download one bundle file; verify size and SHA-256 when the bundle lists them.
ae_download_file <- function(cfg, token, task_id, f, dest) {
  part <- paste0(dest, ".part")
  for (i in 1:3) {
    ok <- tryCatch({
      ae_req(cfg, paste0("bundle/", task_id, "/", f$file_id), token) |>
        httr2::req_perform(path = part)
      TRUE
    }, interrupt = function(e) stop("interrupted by user", call. = FALSE),
       error = function(e) { log_warn("download ", f$file_name, ": ", conditionMessage(e)); FALSE })
    if (ok) {
      size_ok <- is.na(f$file_size) || file.size(part) == f$file_size
      # unclass(): openssl keeps a "hash" class on the hex string, which breaks identical()
      sha_ok <- is.na(f$sha256) ||
        identical(tolower(unclass(as.character(openssl::sha256(file(part))))), tolower(f$sha256))
      if (size_ok && sha_ok) { file.rename(part, dest); return(invisible(TRUE)) }
      log_warn(f$file_name, ": size/checksum mismatch, retrying")
    }
    unlink(part); Sys.sleep(10 * i)
  }
  stop("failed to download ", f$file_name, " after 3 attempts")
}
