#' A dataset about daily average river discharges (in m^3 / s) for the Elbe river network at different measurement stations in Germany
#' 
#' @description Measurements and geographical information about daily average river discharges in (m^3/s) at 13 measurement stations from the Elbe river network from 31.12.1988 to 30.12.2010 
#' for the train data and from 01.01.2010 to 31.12.2020 for the test data.
#' @format A named list containing differnent data files
#' \describe{
#' \item{train}{A list containing the date of the measurement and measurements of the raw discharge data as \code{data.frame} at the 13 stations, 
#' and a \code{data.frame} containing the maximal discharge between the date "from" and "to". 
#' The blockmax dataset only considers the maximal value for the summer months June to September to reduce seasonal trends and temporal dependence.  }
#' \item{test}{Same structure as the two train data.frame objects, but only contains data from 01.01.2011 to 31.12.2020.}
#' \item{info}{A data.frame  object containing the station name, approximate latitude and longitude of the measurement station, the river measured and the next downstream station}
#' }
#' @name elbe
#' @usage data(elbe)
#' @source Datenportal der FGG Elbe <https://www.elbe-datenportal.de>
NULL
