writeNanoStringRccSet <-
function(x, dir = getwd())
{
  # Validate input
  stopifnot(is(x, "NanoStringRccSet"))
  validObject(x)

  # Create directory, if necessary
  if (!dir.exists(dir))
    dir.create(dir)

  # Extract common information
  geneRlf <- annotation(x)
  features <- pData(featureData(x))[, c("CodeClass", "GeneName", "Accession")]
  names(features)[2L] <- "Name"

  # Create output templates
  header1.7 <- "<Header>\nFileVersion,%s\nSoftwareVersion,%s\n</Header>\n"
  header2.0 <- "<Header>\nFileVersion,%s\nSoftwareVersion,%s\nSystemType,%s\n</Header>\n"
  sampleAttr1.7 <-
    paste0("<Sample_Attributes>\nID,%s\nOwner,%s\nComments,%s\nDate,%s\n",
           "GeneRLF,%s\nSystemAPF,%s\n</Sample_Attributes>\n")
  sampleAttr2.0 <-
    paste0("<Sample_Attributes>\nID,%s\nOwner,%s\nComments,%s\nDate,%s\n",
           "GeneRLF,%s\nSystemAPF,%s\nAssayType,%s\n</Sample_Attributes>\n")
  laneAttr <-
    paste0("<Lane_Attributes>\nID,%d\nFovCount,%d\nFovCounted,%d\n",
           "ScannerID,%s\nStagePosition,%d\nBindingDensity,%.2f\n",
           "CartridgeID,%s\nCartridgeBarcode,%s\n</Lane_Attributes>\n")

  # Loop over samples
  for (i in seq_len(dim(x)[["Samples"]])) {
    # Select row data
    protocolRow <- pData(protocolData(x))[i, ]

    # Open file connection
    fname <- file.path(dir, sampleNames(x)[i])
    if (file.exists(fname)) {
      file.remove(fname)
    }
    con <- file(fname, open = "a")

    # Write Header
    if (protocolRow[["FileVersion"]] < numeric_version("2.0"))
      writeLines(sprintf(header1.7,
                         protocolRow[["FileVersion"]],
                         protocolRow[["SoftwareVersion"]]),
                 con)
    else
      writeLines(sprintf(header2.0,
                         protocolRow[["FileVersion"]],
                         protocolRow[["SoftwareVersion"]],
                         protocolRow[["SystemType"]]),
                 con)

    # Write Sample_Attributes
    if (protocolRow[["FileVersion"]] < numeric_version("2.0"))
      writeLines(sprintf(sampleAttr1.7,
                         protocolRow[["SampleID"]],
                         protocolRow[["SampleOwner"]],
                         protocolRow[["SampleComments"]],
                         format(protocolRow[["SampleDate"]], "%Y%m%d"),
                         geneRlf,
                         protocolRow[["SystemAPF"]]),
                 con)
    else
      writeLines(sprintf(sampleAttr2.0,
                         protocolRow[["SampleID"]],
                         protocolRow[["SampleOwner"]],
                         protocolRow[["SampleComments"]],
                         format(protocolRow[["SampleDate"]], "%Y%m%d"),
                         geneRlf,
                         protocolRow[["SystemAPF"]],
                         protocolRow[["AssayType"]]),
                 con)

    # Write Lane_Attributes
    writeLines(sprintf(laneAttr,
                       protocolRow[["LaneID"]],
                       protocolRow[["FovCount"]],
                       protocolRow[["FovCounted"]],
                       protocolRow[["ScannerID"]],
                       protocolRow[["StagePosition"]],
                       protocolRow[["BindingDensity"]],
                       protocolRow[["CartridgeID"]],
                       protocolRow[["CartridgeBarcode"]]),
               con)

    # Write Code_Summary
    writeLines("<Code_Summary>", con)
    write.csv(cbind(features, Count = exprs(x)[, i]), file = con,
              quote = FALSE, row.names = FALSE)
    writeLines("</Code_Summary>\n", con)

    # Write Messages
    writeLines("<Messages>\n</Messages>\n", con)

    # Close file connection
    close(con)
  }

  # Return file paths
  invisible(file.path(dir, sampleNames(x)))
}
