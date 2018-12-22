test_writeNanoStringRccSet_rlf <- function() {
  datadir <- system.file("extdata", "3D_Bio_Example_Data",
                         package = "NanoStringBase")
  files <- dir(datadir, pattern = "SKMEL.*\\.RCC$")
  rcc <-
    readNanoStringRccSet(file.path(datadir, files),
                         file.path(datadir, "3D_SolidTumor_Sig.rlf"))
  writeNanoStringRccSet(rcc, tempdir())
  for (i in seq_along(files)) {
    checkIdentical(readLines(file.path(datadir, files[i])),
                   readLines(file.path(tempdir(), files[i])))
  }
}
