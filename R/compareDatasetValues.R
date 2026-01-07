#' Compare dataset values
#'
#' @description
#' This function is related to `compareDatasets()`, but instead of analysing the
#' motions and outputting graphical results, the aim is to convert the TrackMate
#' XML files into data frames, compile them according to condition and save them
#' in CSV format. Additionally, summary statistics are calculated and saved.
#'
#' Requires TrackMate XML files to be organised into subfolders named according
#' to the condition. The default location for these condition folders is in
#' `Data/` within the working directory. If they are somewhere else, use the
#' `datadir` argument to supply the path. The code will process all the datasets
#' individually, compile them according to condition and compare across
#' conditions. Outputs are saved to `Output/Plots/` in the working directory.
#'
#' If TrackMate XML files require recalibration, this is possible by placing a
#' csv file into each subfolder. All xml files in that folder whose calibration
#' does not match the calibration csv file will be altered. Ideally, all
#' conditions should have the same scaling, and within a condition they should
#' be similar. The code will run if this is not the case, but beware that these
#' discrepancies are not detected. For example, comparing two datasets in um/s
#' with one in mm/min.
#'
#' @param datadir string path to the folder containing condition subfolders with
#'   TrackMate XML files. Default is "Data"
#'
#' @return csv files saved to `Output/Data/`
#' @export
compareDatasetValues <- function(datadir = "Data") {
  condition <- value <- dataid <- cumulative_distance <- NULL
  track_duration <- mean_intensity <- calibrationDF <- timeRes <- NULL
  frame_duration <- speed <- NULL
  units_vec <- NULL
  summary_params <- NULL

  # loop through condition folders within data folder
  condFolderNames <- list.dirs(path = datadir, recursive = FALSE)
  # break if there were no folders in Data directory
  if (identical(condFolderNames, character(0)) == TRUE) {
    # exit with message
    stop("No condition folders found in Data/ directory.
         Please organise TrackMate XML files into subfolders named
         according to condition.")
  }
  # sort the folders alphabetically
  condFolderNames <- sort(condFolderNames)

  for (i in 1:length(condFolderNames)) {
    condFolderPath <- condFolderNames[i]
    condFolderName <- basename(condFolderPath)
    allTrackMateFiles <- list.files(condFolderPath, pattern = "\\.xml$")
    # skip if there were no XML files in this folder
    if (identical(allTrackMateFiles, character(0)) == TRUE) {
      next
    }
    # sort the files alphabetically
    allTrackMateFiles <- sort(allTrackMateFiles)
    # check to see if a calibration file is present
    calibrationFiles <- list.files(condFolderPath, pattern = "\\.csv$")
    calibrationFile <- paste0(condFolderPath, "/", calibrationFiles[1])
    if (identical(calibrationFiles, character(0)) == TRUE) {
      calibrationXY <- 1
      calibrationT <- 1
      calibrate <- FALSE
    } else {
      # load calibration file and store calibrations
      calibDF <- read.csv(calibrationFile)
      calibrate <- TRUE
    }
    cat(paste0("\n", "Processing ", condFolderName, "\n"))
    # Use lists for accumulation
    bigtm_list <- list()
    bigcalibration_list <- list()

    # Parallelize inner loop over XML files
    process_file <- function(j) {
      fileName <- allTrackMateFiles[j]
      thisFilePath <- paste0(condFolderPath, "/", fileName)
      tmObj <- readTrackMateXML(XMLpath = thisFilePath, slim = TRUE)
      if (is.null(tmObj)) {
        cat(paste0("Skipping ", fileName, " - no data found!\n"))
        return(NULL)
      }
      if (calibrate) {
        calibrationDF <- tmObj[[2]]
        calibrationXY <- calibDF[1, 1] / calibrationDF[1, 1]
        calibrationT <- calibDF[2, 1] / calibrationDF[2, 1]
        calibrationXY <- ifelse(calibrationXY == 0, 1, calibrationXY)
        calibrationT <- ifelse(calibrationT == 0, 1, calibrationT)
        calibrationXY <- ifelse(
          calibrationXY < 1.025 & calibrationXY > 0.975, 1, calibrationXY)
        calibrationT <- ifelse(
          calibrationT < 1.025 & calibrationT > 0.975, 1, calibrationT)
        if (calibrationXY != 1 && calibrationT != 1) {
          tmObj <- correctTrackMateData(dataList = tmObj,
                                        xyscalar = calibrationXY,
                                        tscalar = calibrationT,
                                        xyunit = calibDF[1, 2],
                                        tunit = calibDF[2, 2])
        } else if (calibrationXY != 1 && calibrationT == 1) {
          tmObj <- correctTrackMateData(dataList = tmObj,
                                        xyscalar = calibrationXY,
                                        xyunit = calibDF[1, 2])
        } else if (calibrationXY == 1 && calibrationT != 1) {
          tmObj <- correctTrackMateData(dataList = tmObj,
                                        tscalar = calibrationT,
                                        tunit = calibDF[2, 2])
        } else {
          tmObj <- correctTrackMateData(dataList = tmObj,
                                        xyunit = calibDF[1, 2],
                                        tunit = calibDF[2, 2])
        }
      }
      tmDF <- tmObj[[1]]
      calibrationDF <- tmObj[[2]]
      if (is.null(tmDF)) {
        cat(paste0("Skipping ", fileName, " - no data found!\n"))
        return(NULL)
      }

      units_vec <- calibrationDF$unit[1:2]
      thisdataid <- paste0(condFolderName, "_", as.character(j))
      tmDF$dataid <- thisdataid

      return(list(tmDF = tmDF,
                  calibrationDF = calibrationDF))
    }

    # Use mclapply for parallel processing (on non-Windows)
    results <- if (.Platform$OS.type == "windows") {
      lapply(seq_along(allTrackMateFiles), process_file)
    } else {
      mclapply(seq_along(allTrackMateFiles),
               process_file,
               mc.cores = detectCores() - 2)
    }

    for (res in results) {
      if (is.null(res)) next
      if (!is.null(res$tmDF))
        bigtm_list[[length(bigtm_list) + 1]] <- res$tmDF
      if (!is.null(res$calibrationDF))
        bigcalibration_list[[length(bigcalibration_list) + 1]] <-
          res$calibrationDF
    }

    bigtm <- if (length(bigtm_list) > 0)
      do.call(rbind, bigtm_list) else data.frame()

    # Extract units_vec from first valid result
    if (is.null(units_vec)) {
      for (res in results) {
        if (!is.null(res) && !is.null(res$calibrationDF)) {
          units_vec <- res$calibrationDF$unit[1:2]
          break
        }
      }
    }
    # package into a data frame that can be used by makeSummaryReport() to get the units
    dummyCalibrationDF <- data.frame(
      value = c(1, 1),
      unit = units_vec
    )

    bigtmObj <- list(bigtm, dummyCalibrationDF)

    # save data as csv
    destinationDir <- paste0("Output/Data/", condFolderName)
    setupOutputPath(destinationDir)
    # save each dataset-level data
    write.csv(bigtm, paste0(destinationDir, "/allTM.csv"), row.names = FALSE)
  }

  # at this point each condition folder has been processed and saved
  cat("\nAll conditions processed. Summary data saved to Output/Data/.\n")
  # now we will read all of that data back in and compile it
  allTMlist <- list.files("Output/Data",
                          pattern = "allTM.csv$",
                          recursive = TRUE, full.names = TRUE)
  df <- do.call(rbind, lapply(allTMlist, function(x) {
    temp <- read.csv(x)
    temp$condition <- basename(dirname(x))
    temp
  }))
  # calculate summary statistics
  summary_df <- df %>%
    group_by(condition, dataid, trace) %>%
    summarise(frame_duration = n(),
              mean_intensity = mean(mean_intensity, na.rm = TRUE),
              speed = max(cumulative_distance, na.rm = TRUE) /
                max(track_duration, na.rm = TRUE)) %>%
    ungroup()

  write.csv(summary_df, "Output/Data/per_track_summary.csv", row.names = FALSE)

  super_summary_df <- summary_df %>%
    group_by(condition, dataid) %>%
    summarise(tracks = n(),
              avg_frame_duration = mean(frame_duration, na.rm = TRUE),
              avg_mean_intensity = mean(mean_intensity, na.rm = TRUE),
              avg_speed = mean(speed, na.rm = TRUE)) %>%
    ungroup()
  write.csv(super_summary_df, "Output/Data/per_cell_summary.csv", row.names = FALSE)
}
