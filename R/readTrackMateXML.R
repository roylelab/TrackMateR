#' Read TrackMate XML output files.
#'
#' Produces a data frame of all spots from filtered tracks, ordered by track number.
#' A warning is generated if the scaling is in pixels rather than real units.
#'
#' @param XMLpath path to the xml file
#' @param slim if TRUE, only the minimum data required to process the tracks is returned
#' @return list of two data frames
#' @examples
#' xmlPath <- system.file("extdata", "ExampleTrackMateData.xml", package="TrackMateR")
#' tmObj <- readTrackMateXML(XMLpath = xmlPath)
#' # get the track data in a data frame
#' tmDF <-  tmObj[[1]]
#' # get the calibration data in a data frame
#' calibrationDF <- tmObj[[2]]
#' @export


readTrackMateXML <- function(XMLpath, slim = FALSE) {
  # Requires xml2 package
  if (!file.exists(XMLpath)) {
    cat("XML file does not exist: ", XMLpath, "\n")
    return(NULL)
  }

  e <- read_xml(XMLpath)
  track_nodes <- xml_find_all(e, ".//Track")
  if (length(track_nodes) == 0) {
    cat("No tracks found in XML file\n")
    return(NULL)
  }
  filtered_nodes <- xml_find_all(e, ".//TrackID")
  spot_nodes <- xml_find_all(e, ".//AllSpots//SpotsInFrame//Spot")
  feature_nodes <- xml_find_all(e, ".//FeatureDeclarations//SpotFeatures//Feature")
  attrName <- c("name", xml_attr(feature_nodes, "feature"))

  # Units and calibration
  model_node <- xml_find_first(e, ".//Model")
  unitVec <- c(xml_attr(model_node, "spatialunits"), xml_attr(model_node, "timeunits"), "widthpixels", "heightpixels", "ntraces", "maxframes")
  image_node <- xml_find_first(e, ".//ImageData")
  valueVec <- c(xml_attr(image_node, "pixelwidth"), xml_attr(image_node, "timeinterval"), xml_attr(image_node, "width"), xml_attr(image_node, "height"))
  calibrationDF <- data.frame(value = as.numeric(c(valueVec, 0, 0)), unit = unitVec)
  calibrationDF[3:4, 1] <- (calibrationDF[3:4, 1] - 1) * calibrationDF[1, 1]
  calibrationDF[2, 2] <- ifelse(calibrationDF[2, 2] == "sec", "s", calibrationDF[2, 2])
  cat("Units are: ", calibrationDF[1, 1], calibrationDF[1, 2], "and", calibrationDF[2, 1], calibrationDF[2, 2], "\n")
  if (unitVec[1] == "pixel") {
    cat("Spatial units are in pixels - consider transforming to real units\n")
  }
  detector_node <- xml_find_first(e, ".//DetectorSettings")
  target_channel <- xml_attr(detector_node, "TARGET_CHANNEL")
  calibrationDF[7, 1] <- as.numeric(target_channel)
  calibrationDF[7, 2] <- "channel"

  if (slim) {
    slimAttr <- c("name", "POSITION_X", "POSITION_Y", "POSITION_Z", "POSITION_T", "FRAME", "MEAN_INTENSITY", paste0("MEAN_INTENSITY_CH", target_channel))
    attrName <- attrName[attrName %in% slimAttr]
  }

  # Extract spot attributes efficiently
  spot_data <- lapply(attrName, function(attr) xml_attr(spot_nodes, attr))
  dtf <- as.data.frame(spot_data, stringsAsFactors = FALSE)
  for (i in 2:length(attrName)) {
    suppressWarnings(dtf[, i] <- as.numeric(dtf[, i]))
  }
  headerNames <- tolower(attrName)
  if (slim) {
    headerNames <- gsub("_ch\\d$", "", headerNames, ignore.case = TRUE)
  }
  headerNames <- gsub("^position\\w", "", headerNames, ignore.case = TRUE)
  names(dtf) <- headerNames

  cat("Matching track data...\n")
  # Preallocate list for track info
  IDtrace_list <- vector("list", length(track_nodes))
  for (i in seq_along(track_nodes)) {
    tr_node <- track_nodes[[i]]
    edge_nodes <- xml_find_all(tr_node, ".//Edge")
    source_ids <- xml_attr(edge_nodes, "SPOT_SOURCE_ID")
    target_ids <- xml_attr(edge_nodes, "SPOT_TARGET_ID")
    IDvec <- unique(c(source_ids, target_ids))
    trace_id <- xml_attr(tr_node, "TRACK_ID")
    traceDF <- data.frame(name = paste0("ID", IDvec), trace = trace_id, stringsAsFactors = FALSE)
    dispVec <- xml_attr(edge_nodes, "DISPLACEMENT")
    speedVec <- xml_attr(edge_nodes, "SPEED")
    # handle missing analyzers
    if (is.null(dispVec)) dispVec <- rep(0, length(target_ids))
    if (is.null(speedVec)) speedVec <- rep(0, length(target_ids))
    dataDF <- data.frame(name = paste0("ID", target_ids), displacement = as.numeric(dispVec), speed = as.numeric(speedVec), stringsAsFactors = FALSE)
    allDF <- merge(traceDF, dataDF, by = "name", all.x = TRUE)
    allDF$displacement[is.na(allDF$displacement)] <- 0
    allDF$speed[is.na(allDF$speed)] <- 0
    IDtrace_list[[i]] <- allDF
  }
  IDtrace <- do.call(rbind, IDtrace_list)

  daten <- merge(IDtrace, dtf, by = "name")
  FTvec <- xml_attr(filtered_nodes, "TRACK_ID")
  daten <- daten[daten$trace %in% FTvec, ]
  daten <- daten[order(daten$trace, daten$t), ]

  cat("Calculating distances...\n")
  cumdist <- numeric(nrow(daten))
  dur <- numeric(nrow(daten))
  cumdist[1] <- 0
  dur[1] <- 0
  startdur <- daten$t[1]
  for (i in 2:nrow(daten)) {
    if (daten$trace[i] == daten$trace[i - 1]) {
      cumdist[i] <- cumdist[i - 1] + daten$displacement[i]
    } else {
      cumdist[i] <- 0
      startdur <- daten$t[i]
    }
    dur[i] <- daten$t[i] - startdur
  }
  daten$cumulative_distance <- cumdist
  daten$track_duration <- dur

  # Check for multiple spots per frame
  b <- aggregate(name ~ trace + frame, daten, length)
  if (sum(b$name) > nrow(b)) {
    cat("Warning: Detected multiple spots per frame for one or more tracks.\nTrackMateR will only process single tracks. Subsequent analysis will likely fail!\n")
  }

  # Offset coordinates if needed
  minx <- min(daten$x)
  miny <- min(daten$y)
  if (minx < 0) daten$x <- daten$x - minx
  if (miny < 0) daten$y <- daten$y - miny
  maxx <- max(daten$x)
  maxy <- max(daten$y)
  if (maxx > calibrationDF[3, 1]) calibrationDF[3, 1] <- ceiling(maxx)
  if (maxy > calibrationDF[4, 1]) calibrationDF[4, 1] <- ceiling(maxy)
  calibrationDF[5, 1] <- length(unique(daten$trace))
  calibrationDF[6, 1] <- max(daten$track_duration) / calibrationDF[2, 1]

  dfList <- list(daten, calibrationDF)
  return(dfList)
}
