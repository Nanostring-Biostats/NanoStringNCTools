writeNanoStringRccSet <- function(x, dir = getwd()) {
    stopifnot(is(x, "NanoStringRccSet"))
    validObject(x)
    if (!dir.exists(dir)) 
        dir.create(dir)
    geneRlf <- annotation(x)
    features <- pData(featureData(x))[, c("CodeClass", "GeneName", "Accession")]
    names(features)[2L] <- "Name"
    header1.7 <- "<Header>\nFileVersion,%s\nSoftwareVersion,%s\n</Header>\n"
    header2.0 <- "<Header>\nFileVersion,%s\nSoftwareVersion,%s\nSystemType,%s\n</Header>\n"
    sampleAttr1.7 <- paste0("<Sample_Attributes>\nID,%s\nOwner,%s\nComments,%s\nDate,%s\n", 
        "GeneRLF,%s\nSystemAPF,%s\n</Sample_Attributes>\n")
    sampleAttr2.0 <- paste0("<Sample_Attributes>\nID,%s\nOwner,%s\nComments,%s\nDate,%s\n", 
        "GeneRLF,%s\nSystemAPF,%s\nAssayType,%s\n</Sample_Attributes>\n")
    laneAttr <- paste0("<Lane_Attributes>\nID,%d\nFovCount,%d\nFovCounted,%d\n", "ScannerID,%s\nStagePosition,%d\nBindingDensity,%.2f\n", 
        "CartridgeID,%s\nCartridgeBarcode,%s\n</Lane_Attributes>\n")
    for (i in seq_len(dim(x)[["Samples"]])) {
        protocolRow <- pData(protocolData(x))[i, ]
        fname <- file.path(dir, sampleNames(x)[i])
        if (file.exists(fname)) {
            file.remove(fname)
        }
        con <- file(fname, open = "a")
        if (protocolRow[["FileVersion"]] < numeric_version("2.0")) 
            writeLines(sprintf(header1.7, protocolRow[["FileVersion"]], protocolRow[["SoftwareVersion"]]), 
                con)
        else writeLines(sprintf(header2.0, protocolRow[["FileVersion"]], protocolRow[["SoftwareVersion"]], 
            protocolRow[["SystemType"]]), con)
        if (protocolRow[["FileVersion"]] < numeric_version("2.0")) 
            writeLines(sprintf(sampleAttr1.7, protocolRow[["SampleID"]], protocolRow[["SampleOwner"]], 
                protocolRow[["SampleComments"]], format(protocolRow[["SampleDate"]], "%Y%m%d"), 
                geneRlf, protocolRow[["SystemAPF"]]), con)
        else writeLines(sprintf(sampleAttr2.0, protocolRow[["SampleID"]], protocolRow[["SampleOwner"]], 
            protocolRow[["SampleComments"]], format(protocolRow[["SampleDate"]], "%Y%m%d"), 
            geneRlf, protocolRow[["SystemAPF"]], protocolRow[["AssayType"]]), con)
        writeLines(sprintf(laneAttr, protocolRow[["LaneID"]], protocolRow[["FovCount"]], 
            protocolRow[["FovCounted"]], protocolRow[["ScannerID"]], protocolRow[["StagePosition"]], 
            protocolRow[["BindingDensity"]], protocolRow[["CartridgeID"]], protocolRow[["CartridgeBarcode"]]), 
            con)
        writeLines("<Code_Summary>", con)
        write.csv(cbind(features, Count = exprs(x)[, i]), file = con, quote = FALSE, row.names = FALSE)
        writeLines("</Code_Summary>\n", con)
        writeLines("<Messages>\n</Messages>\n", con)
        close(con)
    }
    invisible(file.path(dir, sampleNames(x)))
}
