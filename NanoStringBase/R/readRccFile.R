readRccFile <-
function(file)
{
  # Read data from Reporter Code Count (RCC) file
  lines <- trimws(readLines(file))

  # Split data by tags
  tags <- names(.rccMetadata[["schema"]])
  output <- sapply(tags, function(tag)
  {
    bounds <- charmatch(sprintf(c("<%s>", "</%s>"), tag), lines)
    if (anyNA(bounds) || bounds[1L] + 1L >= bounds[2L])
      lines[integer(0)]
    else
      lines[(bounds[1L] + 1L):(bounds[2L] - 1L)]
  }, simplify = FALSE)

  # Convert single row attributes to data.table objects
  for (tag in c("Header", "Sample_Attributes", "Lane_Attributes")) {
    while (length(bad <- grep(",", output[[tag]], invert = TRUE)) > 0L) {
      bad <- bad[1L]
      if (bad == 1L)
        stop(sprintf("%s section has malformed first line", tag))
      fixed <- output[[tag]]
      fixed[bad - 1L] <- sprintf("%s %s", fixed[bad - 1L], fixed[bad])
      output[[tag]] <- fixed[- bad]
    }
    output[[tag]] <- strsplit(output[[tag]], split = ",")
    output[[tag]] <-
      structure(lapply(output[[tag]],
                       function(x) if (length(x) == 1L) "" else x[2L]),
                names = lapply(output[[tag]], `[`, 1L),
                class = "data.frame",
                row.names = basename(file))
  }

  # Coerce numeric_version information in Header
  metadata <- .rccMetadata[["schema"]][["Header"]]
  cols <- c("FileVersion", "SoftwareVersion")
  if (!(all(cols %in% colnames(output[["Header"]]))))
    stop("Header section must contain \"FileVersion\" and \"SoftwareVersion\"")
  output[["Header"]][, cols] <- lapply(output[["Header"]][, cols],
                                       numeric_version)

  # Extract FileVersion for internal checks
  fileVersion <- output[["Header"]][1L, "FileVersion"]
  if (!(fileVersion %in% numeric_version(c("1.7", "2.0"))))
    stop("\"FileVersion\" in Header section must be either 1.7 or 2.0")

  # Check single row attributes
  for (section in c("Header", "Sample_Attributes", "Lane_Attributes")) {
    valid <- .validRccSchema(output[[section]], fileVersion, section)
    if (!isTRUE(valid))
      stop(valid)
  }

  # Coerce date data type in Sample_Attributes
  output[["Sample_Attributes"]][["Date"]] <-
    as.Date(output[["Sample_Attributes"]][["Date"]], format = "%Y%m%d")

  # Coerce numeric data in Lane Attributes
  cols <- c("ID", "FovCount", "FovCounted", "StagePosition")
  output[["Lane_Attributes"]][, cols] <-
    lapply(output[["Lane_Attributes"]][, cols], as.integer)
  output[["Lane_Attributes"]][["BindingDensity"]] <-
    as.numeric(output[["Lane_Attributes"]][["BindingDensity"]])

  # Rename ID columns in Sample and Lane Attributes
  for (ptag in c("Sample", "Lane")) {
    tag <- sprintf("%s_Attributes", ptag)
    names(output[[tag]])[names(output[[tag]]) == "ID"] <- sprintf("%sID", ptag)
  }

  # Rename Sample Owner, Comments, and Date columns
  for (col in c("Owner", "Comments", "Date")) {
    names(output[["Sample_Attributes"]])[
      names(output[["Sample_Attributes"]]) == col] <- sprintf("Sample%s", col)
  }

  # For FileVersion 1.7, add SystemType to Header and AssayType to Sample Attrs
  if (fileVersion == numeric_version("1.7")) {
    output[["Header"]][["SystemType"]] <- "Gen2"
    output[["Sample_Attributes"]][["AssayType"]] <- NA_character_
  }

  # Convert Code_Summary to data.frame object
  if (output[["Code_Summary"]][1L] != "CodeClass,Name,Accession,Count")
    stop("Code_Summary section header is not \"CodeClass,Name,Accession,Count\"")
  output[["Code_Summary"]][1L] <- "BarcodeClass,GeneName,Accession,Count"
  output[["Code_Summary"]] <- paste(output[["Code_Summary"]], collapse = "\n")
  output[["Code_Summary"]] <-
    read.csv(textConnection(output[["Code_Summary"]]),
             colClasses = c(BarcodeClass = "character", GeneName = "character",
                            Accession = "character", Count = "integer"))
  rownames(output[["Code_Summary"]]) <-
    sprintf("%s_%s_%s",
            output[["Code_Summary"]][["BarcodeClass"]],
            output[["Code_Summary"]][["GeneName"]],
            output[["Code_Summary"]][["Accession"]])

  output
}
