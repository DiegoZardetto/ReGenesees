.onAttach<-function(libname, pkgname){
  path <- read.dcf(file = system.file("DESCRIPTION", package = pkgname, lib.loc = libname))
  packageStartupMessage("\n", appendLF = FALSE)
  packageStartupMessage("\n", appendLF = FALSE)
  packageStartupMessage("--------------------------------------------------------\n", appendLF = FALSE)
  packageStartupMessage("> The ReGenesees package has been successfully loaded. <\n", appendLF = FALSE)
  packageStartupMessage("--------------------------------------------------------\n", appendLF = FALSE)
  packageStartupMessage("\n", appendLF = FALSE)
  packageStartupMessage("\n", appendLF = FALSE)
  write.dcf(path)
  packageStartupMessage("\n", appendLF = FALSE)
}
