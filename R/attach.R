.onAttach <- function(...) {

  tos <- paste0(
    cli::col_green(cli::symbol$info),
    ' ',
    "NIEHS's Terms of Service: ",
    cli::col_blue(cli::style_italic(
      cli::style_hyperlink('<https://niehs.github.io/GeoTox/>', 'https://niehs.github.io/GeoTox/')
    ))
  )

  cite <- paste0(
    cli::col_green(cli::symbol$info),
    ' ',
    'Please cite ', cli::col_blue('GeoToxTimeDependent'), ' if you use it! Use `citation(\'GeoToxTimeDependent\')` for details.'
  )

  rlang::inform(
    paste0(tos, '\n', cite),
    class = 'packageStartupMessage'
  )

  .getAssayIntoPkgEnv(silent = FALSE)
  bootstrap_GeoToxData()
}

.onLoad <- function(...) {
  .getAssayIntoPkgEnv(silent = TRUE)
}


.pkgenv <- new.env(parent = emptyenv())


.getAssayIntoPkgEnv <- function(silent = TRUE){
  .pkgenv[['assay']] <- ''
  assay <- tryCatch({
    ctxR::get_all_assays(Server = ctxR::bioactivity_api_server)
  },
  error = function(e){
    message('Unable to retrieve updated assay information. Please check your API key and internet connection')
    print(e)
    return(NULL)
  }
  )
  .pkgenv[['assay']] <- assay
}

bootstrap_GeoToxData <- function() {
  set_GeoToxTimeDependent_option(
    'httk_data_loaded' = FALSE
  )
}
