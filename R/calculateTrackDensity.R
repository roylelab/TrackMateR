#' Calculate density of tracks
#'
#' Calculate for each track, using its starting frame, what is the relative density of tracks.
#' We use a search radius to find how many tracks in the starting frame are neighbours of the track.
#' The number of neighbours is normalised to the search circle, so that a track in the corner of the image with 2 neighbours has a density of 8.
#' Code for calculating search area (intersection between search circle and the image border) is taken from
#' https://petrelharp.github.io/circle_rectangle_intersection/circle_rectangle_intersection.html
#'
#' @param dataList list of a data frame (must include at a minimum - trace (track ID), x, y and frame (in real coords)) and a calibration data frame
#' @param radius numeric variable for search radius (in spatial units of the data)
#' @return data frame
#' @examples
#' xmlPath <- system.file("extdata", "ExampleTrackMateData.xml", package="TrackMateR")
#' tmObj <- readTrackMateXML(XMLpath = xmlPath)
#' tmObj <- correctTrackMateData(dataList = tmObj, xyscalar = 0.04)
#' tdDF <- calculateTrackDensity(dataList = tmObj, radius = 2)
#' @export


calculateTrackDensity <- function(dataList, radius = 1) {
  # base R optimized version
  if (!inherits(dataList, "list")) {
    cat("Function requires a list of TrackMate data and calibration data\n")
    return(NULL)
  }
  df <- dataList[[1]]
  calibration <- dataList[[2]]
  traceList <- unique(df$trace)

  # Precompute starting positions for all traces
  first_frame_idx <- match(traceList, df$trace)
  frame0s <- df$frame[first_frame_idx]
  x0s <- df$x[first_frame_idx]
  y0s <- df$y[first_frame_idx]

  # Preallocate result vectors
  neighbours <- numeric(length(traceList))
  fractions <- numeric(length(traceList))

  for (j in seq_along(traceList)) {
    frame0 <- frame0s[j]
    x0 <- x0s[j]
    y0 <- y0s[j]
    # select all tracks in the same starting frame
    idx <- which(df$frame == frame0)
    x_all <- df$x[idx]
    y_all <- df$y[idx]
    # vectorized distance calculation
    dists <- sqrt((x_all - x0)^2 + (y_all - y0)^2)
    neighbours[j] <- sum(dists <= radius, na.rm = TRUE) - 1
    # search area fraction
    fractions[j] <- find_td_area(r = radius, xy = c(x0, y0), a = c(0, calibration[3,1]), b = c(0, calibration[4,1])) / (pi * radius^2)
  }

  dfall <- data.frame(trace = traceList,
                      neighbours = neighbours,
                      fraction = fractions)
  dfall$density <- dfall$neighbours / dfall$fraction
  return(dfall)
}

#' Find distance between xy coordinate and a series of other xy coordinates
#'
#' @param xx x coord of point for comparison
#' @param yy y coord of point for comparison
#' @param df data frame containing x and y columns for other points
#' @return numeric vector of distances
#' @keywords internal
find_distances <- function(xx, yy, df) {
  df$x <- (df$x - xx)^2
  df$y <- (df$y - yy)^2
  out <- sqrt(df$x + df$y)

  return(out)
}

#' Track Density - Find area A1
#'
#' @param x value
#' @param r radius
#' @return numeric variable
#' @keywords internal
find_td_A1 <- function (x, r) {
  out <- 0
  if (x < r) {
    out <- r^2 * acos(x/r) - x * sqrt(r^2 - x^2)
  }
  return(out)
}

#' Track Density - Find area A2
#'
#' @param x value
#' @param y value
#' @param r radius
#' @return numeric variable
#' @keywords internal
find_td_A2 <- function (x, y, r) {
  out <- 0
  if (x^2 + y^2 < r^2) {
    out <- (
      (r^2 / 2) * (acos(y/r) + acos(x/r) - pi/2)
      - x * sqrt(r^2 - x^2) / 2
      - y * sqrt(r^2 - y^2) / 2
      + x * y
    )
  }
  return(out)
}

#' Track Density - Find search area
#'
#' Find the area of the intersection of the circle centered at xy with radius r
#' and the radius with vertical sides at a and horizontal sides at b.
#' xy, a, and b must be vectors of length 2, and xy must lie within the rectangle.
#'
#' @param r radius of search circle
#' @param xy numeric vector (length 2)
#' @param a numeric vector (length 2)
#' @param b numeric vector (length 2)
#' @return numeric variable
#' @examples
#' find_td_area(r=2, xy=c(4, 4), a=c(0, 8), b=c(0, 5))
#' @export
#' @keywords internal
find_td_area <- function(r, xy, a, b) {
  stopifnot(length(xy) == 2 && length(a) == 2 && length(b) == 2)
  x1 <- xy[1] - a[1]
  x2 <- a[2] - xy[1]
  y1 <- xy[2] - b[1]
  y2 <- b[2] - xy[2]
  stopifnot(min(x1, x2, y1, y2) >= 0)
  A <- (
    find_td_A1(x1, r)
    + find_td_A1(x2, r)
    + find_td_A1(y1, r)
    + find_td_A1(y2, r)
    - find_td_A2(x1, y1, r)
    - find_td_A2(x1, y2, r)
    - find_td_A2(x2, y1, r)
    - find_td_A2(x2, y2, r)
  )
  stopifnot(A >= 0)
  stopifnot(A <= pi * r^2)

  return(pi * r^2 - A)
}
