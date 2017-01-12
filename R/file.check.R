file.check <- function(file.path) {
    
    if (!(file.access(file.path, 0) == 0)) {
        stop("File '",file.path,"' does not exist.", call. = FALSE)
    }
    
    if (!(file.access(file.path, 4) == 0)) {
        stop("File '",file.path,"' does not have read permission.", call. = FALSE)
    }
}
