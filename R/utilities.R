#' @title Calculate Prevalence Percentage
#' @description Calculate Prevalence Percentage
#' @param dt.vector Data vector
#' @return Prevalence percentage
#' @export
PerctCal <- function(dt.vector) {
  sum(dt.vector!=0)/length(dt.vector)
}


#' @title Calculate Heterogeneity
#' @description Computes heterogeneity using C++ implementation for efficiency
#' @param dt_vector Numeric vector of element values
#' @return Heterogeneity as double
#' @export
#' @examples
#' testdata <- c(1.2, 3.4, 2.1, 5.6)
#' h_val <- HCalC(testdata)
HCalC <- function(dt_vector) {
  if(!is.numeric(dt_vector)) {
    stop("Input must be numeric vector")
  }
  result <- .Call('_qSCMA_HCalC',  dt_vector)
  return(result)
}

#' @title Quantile Calculation Function
#' @description Quantile Calculation Function
#' @param filedir File directory
#' @param num_quan Number of quantiles
#' @param q.cut Quantile cutoff
#' @param heterog Whether to calculate heterogeneity
#' @param prevelen Whether to calculate prevalence
#' @return Reshaped data with quantile features
#' @importFrom dplyr %>%
#' @export
QuanCal <- function(filedir,
                    num_quan,
                    q.cut,
                    heterog,
                    prevelen) {
  tictoc::tic()
  vroom::vroom(filedir,
               progress = FALSE,
               show_col_types = FALSE) -> data

  data %>% names() -> Elements

  batchfunc <- function(Element,
                        data,
                        num_quan,
                        q.cut,
                        heterog,
                        prevelen){
    message("Calculating ",Element)

    data %>%
      dplyr::select(tidyselect::all_of(Element)) %>%
      purrr::set_names("value") %>%
      purrr::as_vector() -> dt.vector

    threshold <- stats::quantile(dt.vector, probs = q.cut/100)

    q.cut.vector <- dt.vector[dt.vector >= threshold]

    q.cut.vector %>%
      stats::quantile(probs = seq(0,1,1/num_quan)) %>%
      as.matrix() %>%
      t() %>% dplyr::as_tibble() %>%
      purrr::set_names(stringr::str_c(Element,
                                      "_Q",
                                      seq(0,1,1/num_quan)*100)) -> QDT

    dt.vector %>%
      mean() %>%
      tibble::tibble(value=.) %>%
      purrr::set_names(stringr::str_c(Element,"_Mean")) %>%
      dplyr::bind_cols(.,
                       QDT) -> QDT


    if(prevelen==T){
      dt.vector %>%
        PerctCal() %>%
        tibble::tibble(value=.) %>%
        purrr::set_names(stringr::str_c(Element,"_P")) %>%
        dplyr::bind_cols(.,
                         QDT) -> QDT
    }

    if(heterog==T){
      dt.vector %>%
        HCalC() %>%
        tibble::tibble(value=.) %>%
        dplyr::mutate(value=ifelse(is.na(value),0,value)) %>%
        purrr::set_names(stringr::str_c(Element,"_H")) %>%
        dplyr::bind_cols(.,
                         QDT) -> QDT
    }

    return(QDT)
  }

  purrr::map_dfc(Elements,
                 ~batchfunc(Element=.,
                            data=data,
                            num_quan=num_quan,
                            q.cut=q.cut,
                            heterog=heterog,
                            prevelen=prevelen)) -> RES

  filedir %>%
    stringr::str_extract("(?<=/)[^/]+(?=\\.)") %>%
    tibble::tibble(filename=.) %>%
    dplyr::bind_cols(.,
                     RES) -> RES
  lubridate::now() %>%
    stringr::str_replace_all(":",".") %>%
    stringr::str_replace_all("-",".") -> NowTime
  tictoc::toc()
  message(filedir," has been reshaped!","\n",NowTime)
  return(RES)
}
