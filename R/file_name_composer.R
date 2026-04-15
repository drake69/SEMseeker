name_composer <- function(...)
{

  details <- c(...)
  # details <- c("prova:zozz","a","scrivere","sto belino","di_nome",".csv")
  fname <- paste0(details, collapse="_")
  fname <- gsub(" ","_", fname)
  fname <- gsub(":","_", fname)
  fname <- gsub("__","_", fname)
  fname <- gsub("_[.]",".", fname)
  fname <- gsub("/_","/", fname)
  fname <- gsub("_/","/", fname)
  fname <- gsub("\\.\\.",".", fname)

  # check if the final character is underscore
  if(fname[length(fname)]=="_")
    fname <- fname[1:(length(fname)-1)]

  return(fname)
}
