#' @title Feature Summary
#' @description near zero variance evaluation for different features
#' @param MT MT object
#' @param features Features to analyze (optional)
#' @param outcome Outcome variable name
#' @param covariate Covariate variables (optional)
#' @param SC Whether to include single-cell data (default: FALSE)
#' @param task Task identifier (optional)
#' @return Updated MT object with feature evaluation results
#' @importFrom dplyr %>%
#' @export
FtEvalu <- function(MT,
                    features = NULL,
                    outcome,
                    covariate = NULL,
                    SC = FALSE,
                    task = NULL) {

  tictoc::tic()

  MT$ReshapedDT -> ReshapedDT

  MT$ReshapedFeatures %>%
    stringr::str_extract(.,"(?<=_)[a-zA-Z]") %>%
    unique() %>% stats::na.omit() %>%
    as.vector() -> FeatureType

  MT$ReshapedFeatures %>%
    tibble::tibble(value=.) -> FT.ALL

  selectfunc <- function(key){
    FT.ALL %>%
      dplyr::filter(stringr::str_detect(value,key)) %>%
      as.matrix() %>% as.vector()
  }

  zerovarfunc <- function(vars,
                          data){
    data %>%
      dplyr::select(tidyselect::all_of(vars)) %>%
      caret::nearZeroVar(saveMetrics = T,
                         names = T) %>%
      dplyr::filter(nzv == "TRUE") %>%
      rownames()  -> vars_deleted

    vars_remain <- setdiff(vars,vars_deleted)

    return(vars_remain)

  }
  eleremainfunc <- function(vars,
                            remain_vars){
    remain_vars %>%
      stringr::str_extract(.,".*(?=_)") %>%
      unique() %>% length() -> n_remain
    vars %>%
      stringr::str_extract(.,".*(?=_)") %>%
      unique() %>% length() -> n_all
    remain <- n_remain/n_all
    return(remain)
  }

  tibble::tibble(FeatureType) %>%
    dplyr::mutate(key=stringr::str_c("_",FeatureType))  %>%
    dplyr::mutate(vars=purrr::map(key,
                                  selectfunc)) %>%
    dplyr::mutate(remain_vars=purrr::map(vars,
                                         ~zerovarfunc(.,
                                                      data=ReshapedDT))) %>%
    dplyr::mutate(remain=purrr::map2_dbl(vars,
                                         remain_vars,
                                         ~eleremainfunc(.x,.y))) %>%
    dplyr::select(-key) -> remaintable

  ## SC -----------------------------------------------------------
  if(SC==T){
    message("Processing single-cell data may require substantial time...")
    MT$RawData$dt_cell -> data.sc
    data.sc %>%
      dplyr::select(-dir.path) %>%
      caret::nearZeroVar(saveMetrics = T,
                         names = T) %>%
      dplyr::filter(nzv == "TRUE") %>%
      rownames()  -> scvars_deleted

    scvars_remain <- setdiff(data.sc %>%
                               dplyr::select(-dir.path) %>%
                               names(),
                             scvars_deleted)
    scremain <- length(scvars_remain)/length(c(scvars_remain,scvars_deleted))

    tibble::tibble(FeatureType="SC",
                   vars=data.sc %>%
                     dplyr::select(-dir.path) %>%
                     names() %>% list(),
           remain_vars=list(scvars_remain),
           remain=scremain) %>%
      dplyr::bind_rows(remaintable,
                       .) -> remaintable
  }



  remaintable %>%
    ggplot2::ggplot(ggplot2::aes(x=FeatureType,y=remain))+
    ggplot2::geom_col(width = 0.75)+
    ggplot2::labs(x="Feature type",
                  y="Non-zero-variance \nelement proportion")+
    ggplot2::theme_classic() -> remain_plot

  MT$FeaturePerformance$table <- remaintable
  MT$FeaturePerformance$plot <- remain_plot

  remaintable %>%
    dplyr::select(-vars) %>%
    dplyr::select(-remain_vars) %>%
    write.csv(.,
              "output/RemainTable.csv")

  lubridate::now() %>%
    stringr::str_replace_all(":",".") %>%
    stringr::str_replace_all("-",".") -> NowTime
  message("Feature summary: Non-near-zero variance\n",
          NowTime)
  print(remaintable)
  print(remain_plot)
  tictoc::toc()
  return(MT)
}

#' @title Feature Identification using Random Forest
#' @description Identify the feature importance by RF model
#' @param MT MT object
#' @param features Features to identify (optional)
#' @param outcome Outcome variable name
#' @param covariate Covariate variables (optional)
#' @param task Task identifier (optional)
#' @return Parameter table for feature identification
#' @importFrom dplyr %>%
#' @export
FtIdent <- function(MT,
                    features = NULL,
                    outcome,
                    covariate = NULL,
                    task = NULL) {
  ## predefine functions for batch process -----------------------------------------------------------
  RF.Boostrap <- function(boostrap,learner.rf,y,x,C,
                          Var.x,data,Path){
    tictoc::tic()
    lubridate::now() %>%
      stringr::str_c("\n",.) -> NowTime
    message(stringr::str_c("Calculating ",y," ~ ",Var.x," for boostrap: ",boostrap,
                           NowTime))
    x <- x
    y <- y
    C <- C

    rows_train <- sample(1:nrow(data),nrow(data)*3/4)
    data %>%
      dplyr::select(tidyselect::all_of(y),
                    tidyselect::all_of(x),
                    tidyselect::all_of(C)) %>%
      .[rows_train,] %>%
      mlr3::as_task_regr(target = y,
                         id = stringr::str_c(y,"~",Var.x)) -> task_train
    data %>%
      dplyr::select(tidyselect::all_of(y),
                    tidyselect::all_of(x),
                    tidyselect::all_of(C)) %>%
      .[-rows_train,] %>%
      mlr3::as_task_regr(target = y,
                         id = stringr::str_c(y,"~",Var.x)) -> task_test

    learner.rf$train(task_train)
    # measure
    measures_train <- list(mlr3::msr("regr.rsq", id = "rsq_train"),
                           mlr3::msr("regr.mae", id = "mae_train"),
                           mlr3::msr("regr.mse", id = "mse_train"),
                           mlr3::msr("regr.rmse", id = "rmse_train"))

    measures_test <- list(mlr3::msr("regr.rsq", id = "rsq_test"),
                          mlr3::msr("regr.mae", id = "mae_test"),
                          mlr3::msr("regr.mse", id = "mse_test"),
                          mlr3::msr("regr.rmse", id = "rmse_test"))
    learner.rf$predict(task_train)$score(measures_train) %>% t() -> msr_train
    learner.rf$predict(task_test)$score(measures_test) %>% t() -> msr_test
    dplyr::bind_cols(Boostrap = boostrap,msr_train,msr_test) -> msr_return

    learner.rf$importance() %>% sort(decreasing=T) -> importance_rf
    importance_rf %>% names() %>% tibble::tibble(vars = .) %>%
      dplyr::mutate(importance = importance_rf) %>%
      dplyr::filter(importance > 0)  %>%
      dplyr::mutate(importance_perct = purrr::map_dbl(importance,
                                                      ~(./sum(importance))),
                    importance_perct_cum = cumsum(importance_perct)) %>%
      purrr::set_names(c("vars",
                         stringr::str_c("importance",boostrap),
                         "importance_perct","importance_perct_cum")) -> importance_table
    importance_table %>%
      writexl::write_xlsx(stringr::str_c(Path,"/importance_",boostrap,".xlsx"))
    tictoc::toc()
    return(msr_return)
  }

  batch.func <- function(Var.y,
                         Var.x,
                         Var.c,
                         dt,
                         task){
    Path <- stringr::str_c("output/FeatureEvaluation_",task,"/",Var.x)
    if(!file.exists(Path)){
      dir.create(Path,recursive = TRUE)}
    ## parameter tuning -----------------------------------------------------------
    lubridate::now() %>%
      stringr::str_c("\n",.) -> NowTime
    message(stringr::str_c("Parameter tunring: ",Var.y," ~ ",Var.x,
                           NowTime))

    dt %>%
      names() %>%
      stringr::str_detect(stringr::str_c(Var.x,"_")) %>%
      dt[,.] %>% names() -> vars
    data_temp <- dt
    data_temp %>%
      dplyr::select(tidyselect::all_of(vars)) %>%
      caret::nearZeroVar(saveMetrics = T,
                         names = T) %>%
      dplyr::filter(nzv == "TRUE") %>%
      rownames()  -> varsdeleted
    x <- setdiff(vars,varsdeleted)
    y <- Var.y
    C <- Var.c

    data_temp %>%
      dplyr::select(tidyselect::all_of(y),
             tidyselect::all_of(x),
             tidyselect::all_of(C)) %>%
      stats::na.omit() -> data_temp
    data_temp %>%
      mlr3::as_task_regr(target = y,
                         id = stringr::str_c(y,"~",Var.x)) -> task_train

    learner.rf = mlr3learners::LearnerRegrRanger$new()
    learner.rf$predict_sets = c("train", "test")

    search_space = paradox::ps(
      max.depth = paradox::p_int(lower = 1, upper = 20),
      mtry.ratio = paradox::p_dbl(lower = 0.2, upper = 0.7),
      num.trees = paradox::p_int(lower = 10, upper = 100)
    )
    terminator = bbotk::trm("evals", n_evals = 500)
    tuner = mlr3tuning::tnr("grid_search")

    at = mlr3tuning::AutoTuner$new(
      learner = learner.rf,
      resampling = mlr3::rsmp("cv", folds = 4),
      measure = mlr3::msr("regr.rsq"),
      search_space = search_space,
      terminator = terminator,
      tuner = tuner
    )
    ddpcr::quiet(at$train(task_train))
    learner.rf$param_set$values <- at$model$tuning_instance$result_learner_param_vals
    learner.rf$param_set$values$importance <- "permutation"


    purrr::map_dfr(1:100,
                   RF.Boostrap,
                   learner.rf = learner.rf,
                   y = y,
                   x = x,
                   C = C,
                   Var.x = Var.x,
                   data = data_temp,
                   Path = Path)  ->  msr_boostrap
    msr_boostrap %>%
      dplyr::mutate(weights1=purrr::map_dbl(rsq_train,
                                            ~ (((./sum(rsq_train))*100))),
                    weights2=purrr::map2_dbl(rsq_train,
                                             rsq_test,
                                             ~ ((((.x+.y)/sum(rsq_train+rsq_test))*100)))) %>%
      writexl::write_xlsx(stringr::str_c(Path,"/msr_boostrap.xlsx"))
    paratable <- tibble::tibble(
      Element = Var.x,
      max.depth = at$model$tuning_instance$result_learner_param_vals$max.depth,
      mtry.ratio = at$model$tuning_instance$result_learner_param_vals$mtry.ratio,
      num.trees = at$model$tuning_instance$result_learner_param_vals$num.trees
    )
    return(paratable)
  }

  ## start processing -----------------------------------------------------------
  if(is.null(task)){
    task <- paste0(
      paste(sample(LETTERS, 5, replace = TRUE),
            collapse = ""),
      paste(sample(0:9, 5, replace = TRUE),
            collapse = "")
    )
  }
  if(!file.exists(stringr::str_c("output/FeatureEvaluation_",task))){
    dir.create(stringr::str_c("output/FeatureEvaluation_",task),recursive = TRUE)}

  MT$ReshapedFeatures %>%
    stringr::str_extract(.,".*(?=_)") %>%
    unique() %>% stats::na.omit() %>%
    as.vector() -> elements

  purrr::map_dfr(elements,
                 ~batch.func(Var.x = .,
                             Var.y=outcome,
                             Var.c=covariate,
                             dt=MT$TransReshapedDT,
                             task = task)) %>%
    writexl::write_xlsx(stringr::str_c(getwd(),"/output/","FeatureEvaluation_",task,"/paratable.xlsx"))


}

