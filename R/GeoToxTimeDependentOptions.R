set_GeoToxTimeDependent_option <- function(...) {

  # if there is no GeoToxData option create the list with the arguments and return
  if (!has_GeoToxData_options()) {
    options('GeoToxData' = list(...))
    return(invisible())
  }

  # otherwise, go through arguments sequentially and add/update
  # them in the list GeoToxData options
  GeoToxData <- getOption('GeoToxData')
  arg_list <- lapply(as.list(match.call())[-1], eval, envir = parent.frame())
  for (k in seq_along(arg_list)) {
    if (names(arg_list)[k] %in% names(GeoToxData)) {
      GeoToxData[names(arg_list)[k]] <- arg_list[k]
    } else {
      GeoToxData <- c(GeoToxData, arg_list[k])
    }
  }

  # set new GeoToxData
  options('GeoToxData' = GeoToxData)

  # return
  invisible()
}


has_GeoToxData_options <- function() {
  !is.null(getOption('GeoToxData'))
}



has_GeoToxData_option <- function(option) {
  if (has_GeoToxData_options()){
    option %in% names(getOption('GeoToxData'))
  } else {
    FALSE
  }
}

httk_data_loaded <- function(){
  if (has_GeoToxData_option('httk_data_loaded')){
    getOption('GeoToxData')$httk_data_loaded
  } else {
    FALSE
  }
}
