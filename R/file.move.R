#' @title Move folders from one location to another
#' @description A small helper function to make file handling
#' in R easier.
#' @param from path from where to copy a folder.
#' @param to path where to copy the corresponding folder specified in \code{from}.
#' @author Hajk-Georg Drost
#' @export

# adapted from: http://stackoverflow.com/questions/10266963/moving-files-between-folders
file.move <- function(from, to) {
  todir <- dirname(to)
  if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
  file.rename(from = from,  to = to)
}