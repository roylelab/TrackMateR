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
  displacement <- NULL
  # Requires xml2 package
  if (!file.exists(XMLpath)) {
    cat("XML file does not exist: ", XMLpath, "\n")
    return(NULL)
  }

  # read the TrackMate XML file
  xml <- read_xml(XMLpath)

  track_nodes <- xml_find_all(xml, ".//Track")
  if (length(track_nodes) == 0) {
    cat("No tracks found in XML file\n")
    return(NULL)
  }

  filtered_nodes <- xml_find_all(xml, ".//TrackID") |>
    xml_attr("TRACK_ID")
  if (length(filtered_nodes) == 0) {
    cat("No filtered tracks found in XML file\n")
    return(NULL)
  }

  # Units and calibration
  model_node <- xml_find_first(xml, ".//Model")
  unitVec <- c(xml_attr(model_node, "spatialunits"),
               xml_attr(model_node, "timeunits"),
               "widthpixels",
               "heightpixels",
               "ntraces",
               "maxframes")
  image_node <- xml_find_first(xml, ".//ImageData")
  valueVec <- c(xml_attr(image_node, "pixelwidth"),
                xml_attr(image_node, "timeinterval"),
                xml_attr(image_node, "width"),
                xml_attr(image_node, "height"))
  calibrationDF <- data.frame(value = as.numeric(c(valueVec, 0, 0)),
                              unit = unitVec)
  calibrationDF[3:4, 1] <- (calibrationDF[3:4, 1] - 1) * calibrationDF[1, 1]
  calibrationDF[2, 2] <- ifelse(calibrationDF[2, 2] == "sec",
                                "s", calibrationDF[2, 2])
  cat("Units are: ", calibrationDF[1, 1], calibrationDF[1, 2], "and",
      calibrationDF[2, 1], calibrationDF[2, 2], "\n")
  if (unitVec[1] == "pixel") {
    cat("Spatial units are in pixels - consider transforming to real units\n")
  }
  # Which channel is the target for tracking by TrackMate?
  detector_node <- xml_find_first(xml, ".//DetectorSettings")
  target_channel <- xml_attr(detector_node, "TARGET_CHANNEL")
  calibrationDF[7, 1] <- as.numeric(target_channel)
  calibrationDF[7, 2] <- "channel"


  spot_nodes <- xml_find_all(xml,
                             ".//AllSpots//SpotsInFrame//Spot")
  feature_nodes <- xml_find_all(xml,
                                ".//FeatureDeclarations//SpotFeatures//Feature")
  attribute_names <- c("name", xml_attr(feature_nodes, "feature"))

  if (slim) {
    slim_attributes <- c("name",
                  "POSITION_X", "POSITION_Y", "POSITION_Z", "POSITION_T",
                  "FRAME", "MEAN_INTENSITY",
                  paste0("MEAN_INTENSITY_CH", target_channel))
    attribute_names <- attribute_names[attribute_names %in% slim_attributes]
  }

  cat("Extracting spot data...\n")

  # Extract spot attribute data efficiently
  spot_data <- extract_spot_attrs(spot_nodes, attribute_names, slim)

  cat("Matching track data...\n")

  # Build spot-to-track mapping (only for filtered tracks)
  spot_track_map <- build_spot_track_map(track_nodes, filtered_nodes)
  spot_track_map$displacement[is.na(spot_track_map$displacement)] <- 0
  spot_track_map$speed[is.na(spot_track_map$speed)] <- 0

  # Join and arrange
  result <- spot_track_map |>
    inner_join(spot_data, by = "name") |>
    select(-"name")
  result <- result[order(as.numeric(result$trace),
                        as.numeric(result$frame)), ]

  # Check for duplicates
  dupe_count <- sum(duplicated(result[, c("trace", "frame")]))
  if (dupe_count > 0) {
    cat("Warning: Detected", dupe_count,
        "duplicate track-frame combinations.\n",
        "TrackMateR will only process single tracks.
        Subsequent analysis will likely fail!\n")
  }

  cat("Calculating distances...\n")

  # Calculate cumulative distance and elapsed time (track duration)
  result <- result %>%
    group_by(trace) %>%
    arrange(frame, .by_group = TRUE) %>%
    mutate(cumulative_distance = cumsum(displacement),
           track_duration = (t - min(t))) %>%
    ungroup()

  # Offset coordinates if necessary
  minx <- min(result$x)
  miny <- min(result$y)
  if (minx < 0) result$x <- result$x - minx
  if (miny < 0) result$y <- result$y - miny
  maxx <- max(result$x)
  maxy <- max(result$y)
  if (maxx > calibrationDF[3, 1]) calibrationDF[3, 1] <- ceiling(maxx)
  if (maxy > calibrationDF[4, 1]) calibrationDF[4, 1] <- ceiling(maxy)
  calibrationDF[5, 1] <- length(unique(result$trace))
  calibrationDF[6, 1] <- max(result$track_duration) / calibrationDF[2, 1]

  dfList <- list(result, calibrationDF)
  return(dfList)
}


#' Extract spot attributes efficiently
#'
#' @param spot_nodes XML nodeset of Spot elements.
#' @param cols list of attributes to extract.
#' @param slim if TRUE, extract only essential columns.
#'
#' @return data frame of spot attributes.
#' @noRd
extract_spot_attrs <- function(spot_nodes, cols, slim) {
  # Pull all attributes at once - much faster than multiple xml_attr calls
  all_attrs <- xml2::xml_attrs(spot_nodes)

  # Extract as matrix then convert
  mat <- vapply(all_attrs, function(x) x[cols], character(length(cols)))
  spots <- as.data.frame(t(mat), stringsAsFactors = FALSE)

  # if a column header is NA, drop it
  na_cols <- which(is.na(names(spots)))
  if (length(na_cols) > 0) {
    spots <- spots[, -na_cols, drop = FALSE]
  }

  # Type conversion
  num_cols <- setdiff(names(spots), "name")
  spots[num_cols] <- lapply(spots[num_cols], as.numeric)
  spots$FRAME <- as.integer(spots$FRAME)

  #
  header_names <- tolower(names(spots))
  if (slim) {
    header_names <- gsub("_ch\\d$", "", header_names, ignore.case = TRUE)
  }
  header_names <- gsub("^position\\w", "", header_names, ignore.case = TRUE)
  names(spots) <- header_names


  return(spots)
}


#' Build a mapping from spot IDs to track IDs
#'
#' @param track_nodes XML nodeset of Track elements.
#' @param filtered_ids Character vector of filtered track IDs to include.
#'
#' @return A data.frame with spot_id and track_id columns.
#' @noRd
build_spot_track_map <- function(track_nodes, filtered_ids) {
  # Pre-filter to only process tracks we care about
  track_ids <- xml2::xml_attr(track_nodes, "TRACK_ID")
  keep <- track_ids %in% filtered_ids
  track_nodes <- track_nodes[keep]
  track_ids <- track_ids[keep]

  # Process each track
  lapply(seq_along(track_nodes), function(i) {
    tr_node <- track_nodes[[i]]
    edge_nodes <- xml2::xml_find_all(tr_node, ".//Edge")
    source_ids <- xml2::xml_attr(edge_nodes, "SPOT_SOURCE_ID")
    target_ids <- xml2::xml_attr(edge_nodes, "SPOT_TARGET_ID")
    IDvec <- unique(c(source_ids, target_ids))

    # Create trace mapping for all spots in track
    trace_data <- data.frame(
      name = paste0("ID", IDvec),
      trace = track_ids[[i]],
      stringsAsFactors = FALSE
    )

    # Extract displacement and speed for target spots
    dispVec <- xml2::xml_attr(edge_nodes, "DISPLACEMENT")
    speedVec <- xml2::xml_attr(edge_nodes, "SPEED")
    if (is.null(dispVec)) dispVec <- rep(0, length(target_ids))
    if (is.null(speedVec)) speedVec <- rep(0, length(target_ids))

    edge_data <- data.frame(
      name = paste0("ID", target_ids),
      displacement = as.numeric(dispVec),
      speed = as.numeric(speedVec),
      stringsAsFactors = FALSE
    )

    # Merge with all.x = TRUE to preserve all spots but add edge data for targets
    merge(trace_data, edge_data, by = "name", all.x = TRUE)
  }) |>
    bind_rows()
}


