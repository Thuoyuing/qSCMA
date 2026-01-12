#' @title Transform Features
#' @description data transformation
#' @param MT MT object
#' @param features Specific features to transform (optional)
#' @param variable Variables to transform (optional)
#' @param logtrans Whether to apply log transformation (default: FALSE)
#' @param scale Whether to apply scaling (default: FALSE)
#' @param factrans Whether to convert to factors (default: FALSE)
#' @return Updated MetaTissue object with transformed data
#' @importFrom dplyr %>%
#' @export
FtTrans <- function(MT,
                    features = NULL,
                    variable = NULL,
                    logtrans = FALSE,
                    scale = FALSE,
                    factrans = FALSE) {
  if(!is.null(features)){
    MT$TransReshapedDT %>%
      dplyr::select(tidyselect::all_of(variable),
                    dplyr::matches(stringr::str_c(stringr::str_c("_",features),
                                                  collapse = "|"))) %>%
      names() -> vars.trans

  } else {
    vars.trans <- variable
  }

  MT$TransReshapedDT -> data

  if(logtrans==T){
    data %>%
      dplyr::mutate(dplyr::across(tidyselect::all_of(vars.trans),
                    ~log(.+0.001))) -> data
  }

  if(scale==T){
    data %>%
      dplyr::mutate(dplyr::across(tidyselect::all_of(vars.trans),
                    ~scale(.,scale = T))) -> data
  }

  if(factrans==T){
    data %>%
      dplyr::mutate(dplyr::across(tidyselect::all_of(vars.trans),
                    ~as.factor(.))) -> data
  }

  data -> MT$TransReshapedDT

  lubridate::now() %>%
    stringr::str_replace_all(":",".") %>%
    stringr::str_replace_all("-",".") -> NowTime

  message("The data has been transformed\n",NowTime)
  return(MT)
}
