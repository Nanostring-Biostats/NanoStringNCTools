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
  header <- "<Header>\nFileVersion,%s\nSoftwareVersion,%s\n</Header>\n"
  sampleAttr <-
    paste0("<Sample_Attributes>\nID,%s\nOwner,%s\nComments,%s\nDate,%s\n",
           "GeneRLF,%s\nSystemAPF,%s\n</Sample_Attributes>\n")
  laneAttr <-
    paste0("<Lane_Attributes>\nID,%d\nFovCount,%d\nFovCounted,%d\n",
           "ScannerID,%s\nStagePosition,%d\nBindingDensity,%.2f\n",
           "CartridgeID,%s\nCartridgeBarcode,%s\n</Lane_Attributes>\n")

  # Loop over samples
  for (i in seq_len(dim(x)[["Samples"]])) {
    # Select row data
    phenoRow <- pData(phenoData(x))[i, ]
    protocolRow <- pData(protocolData(x))[i, ]

    # Open file connection
    fname <- file.path(dir, sampleNames(x)[i])
    if (file.exists(fname)) {
      file.remove(fname)
    }
    con <- file(fname, open = "a")

    # Write Header
    with(protocolRow,
         writeLines(sprintf(header, FileVersion, SoftwareVersion), con))

    # Write Sample_Attributes
    with(phenoRow,
         writeLines(sprintf(sampleAttr, SampleID, Owner, Comments,
                            format(Date, "%Y%m%d"), geneRlf, SystemAPF), con))

    # Write Lane_Attributes
    with(protocolRow,
         writeLines(sprintf(laneAttr, LaneID, FovCount, FovCounted, ScannerID,
                            StagePosition, BindingDensity, CartridgeID,
                            CartridgeBarcode), con))

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
