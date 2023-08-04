.onAttach<-function(libname, pkgname){
  pkg.desc <- read.dcf(file = system.file("DESCRIPTION", package = pkgname, lib.loc = libname))
  conn <- textConnection("RG_Desc_File", "w", local = TRUE)
  sink(conn)
  write.dcf(pkg.desc)
  sink()
  close(conn)
  packageStartupMessage("\n", appendLF = FALSE)
  packageStartupMessage("\n", appendLF = FALSE)
  packageStartupMessage("--------------------------------------------------------\n", appendLF = FALSE)
  packageStartupMessage("> The ReGenesees package has been successfully loaded. <\n", appendLF = FALSE)
  packageStartupMessage("--------------------------------------------------------\n", appendLF = FALSE)
  packageStartupMessage("\n", appendLF = FALSE)
  packageStartupMessage("\n", appendLF = FALSE)
  for (line in RG_Desc_File) {
         packageStartupMessage(line)
    }
  packageStartupMessage("\n", appendLF = FALSE)
}
