### R function to load libraries

## Uses require to check if library is available and then attempts to install 
# package if unavailable

load.library <- function(package) {
  # Load speciifed library
  #
  # Args:
  #   package: character string specifying name of package
  
  # Coerce string into a name
  package <- as.character(substitute(package))
  
  # Load the required package
  if(require(package, character.only = TRUE)) {    
    print(paste(package, "is loaded correctly"))
  } else {
    
    print(paste("trying to install",package))
    install.packages(package)
    
    if(require(package, character.only = TRUE)) {
      print(paste(package, "installed and loaded"))
    } else {
      stop(paste("could not install",package))
    }
  }
  
}
