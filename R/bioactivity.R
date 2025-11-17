

#' Retrieve assay data
#'
#' @param assays A vector of assay names.
#' @param chemids A vector of DTXSID chemical identifiers
#' @param hitc A value between 0 and 1 indicating the cutoff level for assay
#' inclusion
#'
#' @returns A processed tibble with columns `endp` for endpoint, `casn` for
#' CASRN, `dtxsid` for chemical DTXSID, `chnm` for chemical name,
#' `logc` for dose value, and `resp` for assay response.
#' @export
#'
#' @examplesIf FALSE
#' # Retrieve information for `TOX21_DT40` assay for chemicals with DTXSIDs
#' # `DTXSID5020865` and `DTXSID2045078`.
#' g <- get_active_assays(assays =  c('TOX21_DT40'),
#' chemids = c('DTXSID5020865', 'DTXSID2045078'))
#' g
#'
get_active_assays <- function(assays = NULL,
                              chemids = NULL,
                              hitc = 0.9){

  if (is.null(assays) & is.null(chemids)) {
    stop("Must provide at least one 'assays' or 'chemids'")
  }

  hitc_ <- hitc

  # Gather list of assays
  assay_list <- assay_data()
  if (is.null(assay_list)){
    assay_list <- ctxR::get_all_assays()
  }
  # print('Assay list')

  # Select assays given by user input and return AEID assay identifiers
  aeids <- assay_list |> dplyr::filter(.data$assayComponentEndpointName %in% assays) |>
    dplyr::select(.data$aeid, endp = .data$assayName)
  # print('AEIDs')
  # print(aeids)
  # print(names(aeids))

  # Retrieve assay data based on input chemical list
  assay_data <- data.table::rbindlist(ctxR::get_bioactivity_details_batch(DTXSID = chemids))
  # print('Assay data')
  # print(dim(assay_data))
  # print(sum(assay_data$hitc >= hitc_))

  assay_cols <- c('casn', 'dtxsid', 'aeid', 'conc', 'resp', 'spid')
  # Filter assay data to relevant assays by AEID values and hitc values
  assay_data <- assay_data |> dplyr::filter(.data$aeid %in% aeids$aeid) |>
    dplyr::filter(.data$hitc >= hitc_) |>
    dplyr::select(tidyselect::any_of(assay_cols)) |>
    tidyr::unnest(cols = c('conc', 'resp'))
  # print('Filter by aeid')
  # print(dim(assay_data))
  # print(sum(assay_data$hitc >= hitc_))
  #print(assay_data)

  # Join assay data
  assay_data <- assay_data |> dplyr::left_join(aeids, by = dplyr::join_by('aeid'))
  # print('Assaydata')
  # print(names(assay_data))

  # Retrieve chemical details
  chemicals <- ctxR::get_chemical_details_batch(DTXSID = chemids) |>
    dplyr::select(c(.data$dtxsid, chnm = .data$preferredName))

  # Join assay data and chemical data
  assay_data <- assay_data |>
    dplyr::inner_join(chemicals, by = dplyr::join_by('dtxsid')) |> as.data.frame()
  # print(dim(assay_data))


  output_cols <- c('endp', 'casn', 'dtxsid', 'chnm', 'spid')
  assay_data |> tidyr::unnest(cols = c('conc', 'resp')) |>
    dplyr::mutate(c = conc) |> dplyr::mutate(logc = log10(conc)) |> dplyr::select(-.data$aeid) |>
    dplyr::relocate(output_cols)

  # # Process concentration and response data
  # assay_data |> dplyr::rowwise(function(t) {
  # tibble(casn = t['casn'],
  #        dtxsid = t['dtxsid'],
  #        aeid = t['aeid'],
  #        logc = 10^unlist(t['logc']),
  #        resp = unlist(t['resp']),
  #        preferredName = t['preferredName']
  # )
  # }
  # )
}


#' Retrieve assay data from environmental variables
#'
#' @returns data.frame of assay data.
#' @export
#'
#' @examplesIf FALSE
#' # Examine assay data
#' head(b)
assay_data <- function() {
  assays <- .pkgenv[["assay"]]
  return(assays)
}
