ws.wrap.path <- function(paths){
  
  paths <- stringr::str_split(paths,.Platform$file.sep)
  
  ws.paths <- unlist(lapply(paths, function(x) {
    
    stringr::str_c(unlist(sapply(x, function(word) {
      if(stringr::str_detect(word," ")) 
        return(stringr::str_replace(word,word,paste0("'",word,"'"))) else 
          return(word)})),collapse = .Platform$file.sep)
    
  }))
 
  return (ws.paths) 
}