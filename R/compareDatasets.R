#' Compare datasets
#'
#' Requires TrackMate XML files to be organised into subfolders named according to the condition.
#' If these condition folders are in `Data/` within the working directory, the code will run automatically.
#' Since there is no easy cross-platform way for the user to interactively pick a directory, the organisation of files in `Data/` is a requirement.
#' The `Data/` folder can be elsewhere on your computer, just change the `wd` prior to running `compareDatasets()` and the routine will run.
#' The code will process all the datasets individually, compile them according to condition and compare across conditions.
#' Outputs are saved to `Output/Plots/` in the working directory.
#'
#' If TrackMate XML files require recalibration, this is possible by placing a csv file into each subfolder.
#' All xml files in that folder whose calibration does not match the calibration csv file will be altered.
#' Ideally, all conditions should have the same scaling, and within a condition they should be similar.
#' The code will run if this is not the case, but beware that these discrepancies are not detected.
#' For example, comparing two datasets in um/s with one in mm/min.
#'
#' @param ... pass additional parameters to modify the defaults (N, short, deltaT, mode, nPop, init, timeRes, breaks, radius)
#'
#' @return multiple pdf reports
#' @export
compareDatasets <- function(...) {

  condition <- value <- dataid <- cumulative_distance <- track_duration <- mean_intensity <- calibrationDF <- timeRes <-NULL
  megamsd <- megaalpha <- megadee <- megatd <- megaspeed <- megafd <- NULL
  units_vec <- NULL
  summary_params <- NULL

  if(!dir.exists("Data")) {
    # there is no cross-platform way to safely choose directory
    cat("Please organise your XML files in a folder called Data in the working directory\n")
    return(-1)
  } else {
    datadir <- "Data"
  }

  # ellipsis processing
  l <- NULL
  l <- list(...)
  l <- processEllipsis(l)

  # loop through condition folders within data folder
  condFolderNames <- list.dirs(path = datadir, recursive = FALSE)
  # break if there were no folders in Data directory
  if(identical(condFolderNames, character(0)) == TRUE) {
    cat("Please organise your XML files in a subfolders within Data in the
        working directory\n")
    return(-1)
  }
  # sort the folders alphabetically
  condFolderNames <- sort(condFolderNames)

  for(i in 1:length(condFolderNames)) {
    condFolderPath <- condFolderNames[i]
    condFolderName <- basename(condFolderPath)
    allTrackMateFiles <- list.files(condFolderPath, pattern = "*.xml")
    # skip if there were no XML files in this folder
    if(identical(allTrackMateFiles, character(0)) == TRUE) {
      next
    }
    # sort the files alphabetically
    allTrackMateFiles <- sort(allTrackMateFiles)
    # check to see if a calibration file is present
    calibrationFiles <- list.files(condFolderPath, pattern = "*.csv")
    calibrationFile <- paste0(condFolderPath,"/",calibrationFiles[1])
    if(identical(calibrationFiles, character(0)) == TRUE) {
      calibrationXY <- 1
      calibrationT <- 1
      calibrate <- FALSE
    } else {
      # load calibration file and store calibrations
      calibDF <- read.csv(calibrationFile)
      calibrate <- TRUE
    }
    cat(paste0("\n","Processing ",condFolderName,"\n"))
    # Use lists for accumulation
    bigtm_list <- list()
    bigcalibration_list <- list()
    bigmsd_list <- list()
    bigalpha_list <- list()
    bigdee_list <- list()
    bigjd_list <- list()
    bigjdParams_list <- list()
    bigtd_list <- list()
    bigfd_list <- list()
    bigreport_list <- list()

    # Parallelize inner loop over XML files
    process_file <- function(j) {
      fileName <- allTrackMateFiles[j]
      thisFilePath <- paste0(condFolderPath, "/", fileName)
      tmObj <- readTrackMateXML(XMLpath = thisFilePath, slim = TRUE)
      if(is.null(tmObj)) {
        cat(paste0("Skipping ",fileName," - no data found!\n"))
        return(NULL)
      }
      if(calibrate) {
        calibrationDF <- tmObj[[2]]
        calibrationXY <- calibDF[1,1] / calibrationDF[1,1]
        calibrationT <- calibDF[2,1] / calibrationDF[2,1]
        calibrationXY <- ifelse(calibrationXY == 0, 1, calibrationXY)
        calibrationT <- ifelse(calibrationT == 0, 1, calibrationT)
        calibrationXY <- ifelse(calibrationXY < 1.025 & calibrationXY > 0.975, 1, calibrationXY)
        calibrationT <- ifelse(calibrationT < 1.025 & calibrationT > 0.975, 1, calibrationT)
  if(calibrationXY != 1 && calibrationT != 1) {
          tmObj <- correctTrackMateData(dataList = tmObj, xyscalar = calibrationXY, tscalar = calibrationT, xyunit = calibDF[1,2], tunit = calibDF[2,2])
  } else if(calibrationXY != 1 && calibrationT == 1) {
          tmObj <- correctTrackMateData(dataList = tmObj, xyscalar = calibrationXY, xyunit = calibDF[1,2])
  } else if(calibrationXY == 1 && calibrationT != 1) {
          tmObj <- correctTrackMateData(dataList = tmObj, tscalar = calibrationT, tunit = calibDF[2,2])
        } else {
          tmObj <- correctTrackMateData(dataList = tmObj, xyunit = calibDF[1,2], tunit = calibDF[2,2])
        }
      }
      tmDF <- tmObj[[1]]
      calibrationDF <- tmObj[[2]]
      if(is.null(tmDF)) {
        cat(paste0("Skipping ",fileName," - no data found!\n"))
        return(NULL)
      }
      if(calibrationDF[5,1] < 3 & calibrationDF[6,1] < 10) {
        cat(paste0("Skipping ",fileName," as it has less than 3 tracks and longest track has less than 10 frames\n"))
        return(NULL)
      }
  units_vec <- calibrationDF$unit[1:2]
      thisdataid <- paste0(condFolderName,"_",as.character(j))
      tmDF$dataid <- thisdataid
      msdObj <- calculateMSD(tmDF, N = 3, short = 8)
      msdDF <- msdObj[[1]]
      alphaDF <- msdObj[[2]]
      deeDF <- msdObj[[3]]
  if(!is.null(msdDF) || !is.null(alphaDF) || !is.null(deeDF)) {
        msdDF$dataid <- thisdataid
        alphaDF$dataid <- thisdataid
        deeDF$dataid <- thisdataid
      }
      deltaT <- 1
      jdObj <- calculateJD(dataList = tmObj, deltaT = l$deltaT, nPop = l$nPop, mode = l$mode, init = l$init, timeRes = l$timeRes, breaks = l$breaks)
      jdDF <- jdObj[[1]]
      if(is.null(jdDF)) {
        jdObj <- NULL
      } else {
        jdDF$dataid <- thisdataid
        timeRes <- jdObj[[2]]
        jdObj <- list(jdDF,timeRes)
      }
      tdDF <- calculateTrackDensity(dataList = tmObj, radius = l$radius)
      if(!is.null(tdDF)) {
        tdDF$dataid <- thisdataid
      }
      fdDF <- calculateFD(dataList = tmObj)
      if(!is.null(fdDF)) {
        fdDF$dataid <- thisdataid
      }
      fileNameBase <- tools::file_path_sans_ext(basename(thisFilePath))
      both <- makeSummaryReport(tmList = tmObj, msdList = msdObj, jumpList = jdObj, tddf = tdDF, fddf = fdDF,
                                titleStr = condFolderName, subStr = fileNameBase, auto = TRUE, summary = FALSE,
                                msdplot = l$msdplot)
      p <- both[[1]]
      destinationDir <- paste0("Output/Plots/", condFolderName)
      setupOutputPath(destinationDir)
      filePath <- paste0(destinationDir, "/report_",as.character(j),".pdf")
      ggsave(filePath, plot = p, width = 25, height = 19, units = "cm")
      df_report <- both[[2]]
      df_report$condition <- condFolderName
      df_report$dataid <- thisdataid
      return(list(tmDF = tmDF, calibrationDF = calibrationDF, msdDF = msdDF, alphaDF = alphaDF, deeDF = deeDF, jdDF = jdDF, jdParams = timeRes, tdDF = tdDF, fdDF = fdDF, df_report = df_report))
    }

    # Use mclapply for parallel processing (on non-Windows)
    results <- if (.Platform$OS.type == "windows") {
      lapply(seq_along(allTrackMateFiles), process_file)
    } else {
      mclapply(seq_along(allTrackMateFiles), process_file, mc.cores = detectCores() - 2)
    }

    for (res in results) {
      if (is.null(res)) next
      if (!is.null(res$tmDF)) bigtm_list[[length(bigtm_list)+1]] <- res$tmDF
      if (!is.null(res$calibrationDF)) bigcalibration_list[[length(bigcalibration_list)+1]] <- res$calibrationDF
      if (!is.null(res$msdDF)) bigmsd_list[[length(bigmsd_list)+1]] <- res$msdDF
      if (!is.null(res$alphaDF)) bigalpha_list[[length(bigalpha_list)+1]] <- res$alphaDF
      if (!is.null(res$deeDF)) bigdee_list[[length(bigdee_list)+1]] <- res$deeDF
      if (!is.null(res$jdDF)) bigjd_list[[length(bigjd_list)+1]] <- res$jdDF
      if (!is.null(res$jdParams)) bigjdParams_list[[length(bigjdParams_list)+1]] <- res$jdParams
      if (!is.null(res$tdDF)) bigtd_list[[length(bigtd_list)+1]] <- res$tdDF
      if (!is.null(res$fdDF)) bigfd_list[[length(bigfd_list)+1]] <- res$fdDF
      if (!is.null(res$df_report)) bigreport_list[[length(bigreport_list)+1]] <- res$df_report
    }

  bigtm <- if (length(bigtm_list) > 0) do.call(rbind, bigtm_list) else data.frame()
  bigmsd <- if (length(bigmsd_list) > 0) do.call(rbind, bigmsd_list) else data.frame()
  bigalpha <- if (length(bigalpha_list) > 0) do.call(rbind, bigalpha_list) else data.frame()
  bigdee <- if (length(bigdee_list) > 0) do.call(rbind, bigdee_list) else data.frame()
  bigjd <- if (length(bigjd_list) > 0) do.call(rbind, bigjd_list) else data.frame()
  bigtd <- if (length(bigtd_list) > 0) do.call(rbind, bigtd_list) else data.frame()
  bigfd <- if (length(bigfd_list) > 0) do.call(rbind, bigfd_list) else data.frame()
  bigreport <- if (length(bigreport_list) > 0) do.call(rbind, bigreport_list) else data.frame()

  # Fix: collect all valid parameter lists from results and use the first one for summary
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
dummyCalibrationDF <- data.frame(value = c(1, 1),
                             unit = units_vec)

# Extract summary_params from first valid result
all_param_lists <- lapply(results, function(res) {
  if (!is.null(res) && !is.null(res$jdParams)) res$jdParams else NULL
})
valid_param_lists <- Filter(Negate(is.null), all_param_lists)
if (is.null(summary_params) || length(summary_params) == 0) {
  summary_params <- if (length(valid_param_lists) > 0) valid_param_lists[[1]] else list()
}
  bigtmObj <- list(bigtm, dummyCalibrationDF)
  bigmsdObj <- list(bigmsd, bigalpha, bigdee)
  bigjdObj <- list(bigjd, summary_params)
  summaryObj <- makeSummaryReport(tmList = bigtmObj, msdList = bigmsdObj, jumpList = bigjdObj, tddf = bigtd, fddf = bigfd,
                  titleStr = condFolderName, subStr = "Summary", auto = TRUE, summary = TRUE,
                  msdplot = l$msdplot)
  p <- summaryObj[[1]]
  destinationDir <- paste0("Output/Plots/", condFolderName)
  filePath <- paste0(destinationDir, "/combined.pdf")
  ggsave(filePath, plot = p, width = 25, height = 19, units = "cm")
    # save data as csv
    destinationDir <- paste0("Output/Data/", condFolderName)
    setupOutputPath(destinationDir)
    # save each dataset-level data
    write.csv(bigtm, paste0(destinationDir, "/allTM.csv"), row.names = FALSE)
    write.csv(bigmsd, paste0(destinationDir, "/allMSD.csv"), row.names = FALSE)
    write.csv(bigjd, paste0(destinationDir, "/allJD.csv"), row.names = FALSE)
    write.csv(bigfd, paste0(destinationDir, "/allFD.csv"), row.names = FALSE)
    # mega data frame of msd averages per dataset; alpha values, track density, speed/duration/distance, intensity by trace/dataid/condition
    msdSummary <- summaryObj[[2]]
    msdSummary$condition <- condFolderName
    bigspeed <- bigtm %>%
      group_by(dataid, trace) %>%
      summarise(cumdist = max(cumulative_distance), cumtime = max(track_duration), intensity = max(mean_intensity))
    bigspeed$speed <- bigspeed$cumdist / bigspeed$cumtime
    bigspeed$condition <- condFolderName
    if (is.null(megamsd)) {
      megamsd <- msdSummary
      megaalpha <- bigalpha
      megadee <- bigdee
      megatd <- bigtd
      megaspeed <- bigspeed
      megafd <- bigfd
      megareport <- bigreport
    } else {
      megamsd <- rbind(megamsd, msdSummary)
      megaalpha <- rbind(megaalpha, bigalpha)
      megadee <- rbind(megadee, bigdee)
      megatd <- rbind(megatd, bigtd)
      megaspeed <- rbind(megaspeed, bigspeed)
      megafd <- rbind(megafd, bigfd)
      megareport <- rbind(megareport, bigreport)
    }
  }

  # save summary data as csv
  destinationDir <- "Output/Data"
  write.csv(megamsd, paste0(destinationDir, "/allMSDCurves.csv"), row.names = FALSE)
  write.csv(megareport, paste0(destinationDir, "/allComparison.csv"), row.names = FALSE)

  # for alpha values, estimator of D, track density, peed/duration/distance, intensity by trace/dataid/condition we must combine into one
  # set the name of dee to estdee
  names(megadee)[names(megadee) == "dee"] <- "estdee"
  megatrace <- Reduce(mergeDataFramesForExport, list(megaalpha, megadee, megatd, megaspeed, megafd))
  write.csv(megatrace, paste0(destinationDir, "/allTraceData.csv"), row.names = FALSE)

  # generate the comparison plots and save
  p <- makeComparison(df = megareport, msddf = megamsd, units = units_vec, msdplot = l$msdplot)
  destinationDir <- "Output/Plots/"
  filePath <- paste0(destinationDir, "/comparison.pdf")
  # if there are many conditions, increase width of plot
  if(length(condFolderNames) < 3) {
    ggsave(filePath, plot = p, width = 19, height = 19, units = "cm")
  } else if (length(condFolderNames) > 6) {
    ggsave(filePath, plot = p, width = 35, height = 19, units = "cm")
  } else {
    ggsave(filePath, plot = p, width = 25, height = 19, units = "cm")
  }
}
