#' @title Initialize qSCMA Analysis Object
#' @description Creates an R6 object for storing qSCMA analysis data and results
#' @param cell_data_path Path to directory containing cell data files
#' @param outcome_data Path to outcome data file
#' @return An R6 object containing data paths and storage for analysis results
#' @importFrom dplyr %>%
#' @export
#' @examples
#' \dontrun{
#' MT <- ReadMT("path/to/cell/data", "path/to/outcome.csv")
#' }
ReadMT <- function(cell_data_path, outcome_data) {

  tictoc::tic()

  # creat output dir
  Path <- "output/"
  if(!dir.exists(Path)) dir.create(Path, recursive = TRUE)

  # define R6 class
  MT <- R6::R6Class(
    "MT",
    public = list(
      Dir.RawData = list(dir.outc = NULL, dir.cell = NULL),
      RawData = list(dt_cell = NULL, dt_outc = NULL),
      ReshapedDT = NULL,
      TransReshapedDT = NULL,
      FeaturePerformance = NULL
    ),
    lock_class = FALSE,
    lock_objects = FALSE
  )

  # 创建实例
  MT_instance <- MT$new()

  # 设置数据路径
  cell_files <- dir(cell_data_path, full.names = TRUE)
  outcome_data <- vroom::vroom(outcome_data, show_col_types = FALSE)

  MT_instance$Dir.RawData$dir.cell <- cell_files
  MT_instance$Dir.RawData$dir.outc <- outcome_data
  MT_instance$RawData$dt_outc <- outcome_data

  # 记录完成时间
  NowTime <- format(Sys.time(), "%Y.%m.%d.%H.%M.%S")
  message("Data preloaded successfully!\n",
          "Outcome data rows: ", nrow(outcome_data), "\n",
          "Cell data files: ", length(cell_files), "\n",
          NowTime)
  tictoc::toc()
  return(MT_instance)
}

#' @title Load Single Cell Data
#' @description Loading single-cell data to MT
#' @param MT MT
#' @return Updated MT object with loaded cell data
#' @importFrom dplyr %>%
#' @export
LoadSCDT <- function(MT) {

  tictoc::tic()
  MT$Dir.RawData$dir.cell -> dir.cell

  bindfunc <- function(files){
    vroom::vroom(files,
                 show_col_types = FALSE) %>%
      dplyr::mutate(dir.path=files) %>%
      dplyr::select(dir.path,
                    tidyselect::everything()) -> Res
    return(Res)
  }

  purrr::map_dfr(dir.cell,
                 ~bindfunc(.)) -> MT$RawData$dt_cell

  lubridate::now() %>%
    stringr::str_replace_all(":",".") %>%
    stringr::str_replace_all("-",".") -> NowTime

  message("Cell data loaded!\n--number of cells: ",
          nrow(MT$RawData$dt_cell))
  message(NowTime)
  tictoc::toc()
  return(MT)
}
