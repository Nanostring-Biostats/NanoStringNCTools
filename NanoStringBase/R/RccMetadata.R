.rccMetadata <-
  list(schema =
         list("Header" =
                data.frame(labelDescription =
                             c("The version of the file",
                               "The version of the software used to create the file",
                               "The generation of instrument used to create the file"),
                           minVersion = numeric_version(c("1.7", "1.7", "2.0")),
                           row.names =
                             c("FileVersion", "SoftwareVersion", "SystemType"),
                           stringsAsFactors = FALSE),
              "Sample_Attributes" =
                data.frame(labelDescription =
                             c("The sample ID",
                               "The owner of the sample",
                               "Comments about the sample",
                               "The date of the sample",
                               "The filename without the \"*.RLF\" extension",
                               "The filename without the \"*.APF\" extension",
                               "The assay type that the user associated with this codeset"),
                           minVersion = numeric_version(c(rep("1.7", 6L), "2.0")),
                           row.names = 
                             c("ID", "Owner", "Comments", "Date", "GeneRLF",
                               "SystemAPF", "AssayType"),
                           stringsAsFactors = FALSE),
              "Lane_Attributes" =
                data.frame(labelDescription =
                             c("The lane ID (1-12)",
                               "The specified FOV count",
                               "The number of FOV counted",
                               "The ID of the scanner used",
                               "The stage position of the lane's cartridge (1-6)",
                               "The density of spots in the lane",
                               "The ID of the cartridge",
                               "The barcode of the cartridge"),
                           minVersion = numeric_version(c(rep("1.7", 8L))),
                           row.names =
                             c("ID", "FovCount", "FovCounted", "ScannerID",
                               "StagePosition", "BindingDensity", "CartridgeID",
                               "CartridgeBarcode"),
                           stringsAsFactors = FALSE),
              "Code_Summary" =
                data.frame(labelDescription =
                             c(NA_character_, NA_character_, NA_character_,
                               NA_character_),
                           minVersion = numeric_version(c(rep("1.7", 4L))),
                           row.names = c("CodeClass", "Name", "Accession", "Count"),
                           stringsAsFactors = FALSE),
              "Messages" = "character")
  )


.rccMetadata[["protocolData"]] <-
  do.call(rbind,
          unname(head(.rccMetadata[["schema"]], 3L)))[, "labelDescription",
                                                      drop = FALSE]
.rccMetadata[["protocolData"]] <-
  .rccMetadata[["protocolData"]][rownames(.rccMetadata[["protocolData"]]) !=
                                   "GeneRLF", , drop = FALSE]
rownames(.rccMetadata[["protocolData"]])[
  rownames(.rccMetadata[["protocolData"]]) == "ID"] <- "SampleID"
rownames(.rccMetadata[["protocolData"]])[
  rownames(.rccMetadata[["protocolData"]]) == "ID1"] <- "LaneID"


.validRccSchema <-
function(x, fileVersion,
         section = c("Header", "Sample_Attributes", "Lane_Attributes",
                     "Code_Summary"))
{
  section <- match.arg(section)
  schema <- .rccMetadata[["schema"]][[section]]
  expectedNames <- row.names(schema)[schema[,"minVersion"] <= fileVersion]
  if (identical(colnames(x), expectedNames))
    TRUE
  else
    sprintf("<%s> section must contain %s", section,
            paste0("\"", expectedNames, "\"", collapse = ", "))
}

.validNonNegativeInteger <- function(x) {
  is.integer(x) && !anyNA(x) && min(x) >= 0L
}

.validNonNegativeNumber <- function(x) {
  is.numeric(x) && !anyNA(x) && min(x) >= 0
}
