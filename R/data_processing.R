#' @title Reshape Data
#' @description Reshape cell-level data to sample level features
#' @param MT MT object
#' @param num_quan Number of quantiles (default: 100)
#' @param q.cut Quantile cutoff threshold (default: 0)
#' @param heterog Whether to calculate heterogeneity (default: TRUE)
#' @param prevelen Whether to calculate prevalence (default: TRUE)
#' @return Updated MT object with reshaped data
#' @importFrom dplyr %>%
#' @export
DtReshape <- function(MT,
                      num_quan = 100,
                      q.cut = 0,
                      heterog = TRUE,
                      prevelen = TRUE) {

  tictoc::tic()

  MT$Dir.RawData$dir.cell

  purrr::map_dfr(MT$Dir.RawData$dir.cell,
          ~QuanCal(filedir = .,
                   num_quan = num_quan,
                   q.cut = q.cut,
                   heterog = heterog,
                   prevelen = prevelen)) -> QData
  QData %>% utils::write.csv(.,
                             "output/ReshapedDT.CSV")
  QData -> MT$ReshapedDT
  QData %>%
    dplyr::select(-filename) %>%
    names() -> MT$ReshapedFeatures
  lubridate::now() %>%
    stringr::str_replace_all(":",".") %>%
    stringr::str_replace_all("-",".") -> NowTime
  message("Data reshaped!\n",
          NowTime)
  tictoc::toc()
  return(MT)
}

#' @title Bind Data
#' @description Bind reshaped single-cell data to outcome data
#' @param MT MT object
#' @return Updated MT object with outcome data bound
#' @importFrom dplyr %>%
#' @export
DtBind <- function(MT) {
  MT$ReshapedDT %>%
    dplyr::mutate(filename=as.character(filename)) -> MT$ReshapedDT
  MT$RawData$dt_outc %>%
    dplyr::mutate(filename=as.character(filename)) -> MT$RawData$dt_outc

  dplyr::inner_join(MT$RawData$dt_outc,
                    MT$ReshapedDT,
                    by="filename") -> MT$ReshapedDT
  MT$TransReshapedDT <- MT$ReshapedDT
  lubridate::now() %>%
    stringr::str_replace_all(":",".") %>%
    stringr::str_replace_all("-",".") -> NowTime
  message("The ReshapedDT dataset was updated with the outcome data bound.\n",NowTime)
  return(MT)
}
