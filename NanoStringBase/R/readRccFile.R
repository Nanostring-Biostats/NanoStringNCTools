readRccFile <-
function(file)
{
  # Read data from Reporter Code Count (RCC) file
  lines <- trimws(readLines(file))

  # Split data by tags
  tags <- c("Header", "Sample_Attributes", "Lane_Attributes", "Code_Summary",
            "Messages")
  output <- sapply(tags, function(tag)
  {
    bounds <- charmatch(sprintf(c("<%s>", "</%s>"), tag), lines)
    if (anyNA(bounds) || bounds[1L] + 1L >= bounds[2L])
      lines[integer(0)]
    else
      lines[(bounds[1L] + 1L):(bounds[2L] - 1L)]
  }, simplify = FALSE)

  # Convert single row attributes to data.table objects
  for (tag in intersect(c("Header", "Sample_Attributes", "Lane_Attributes"),
                        names(output))) {
    output[[tag]] <- strsplit(output[[tag]], split = ",")
    output[[tag]] <-
      structure(lapply(output[[tag]],
                       function(x) if (length(x) == 1L) "" else x[2L]),
                names = lapply(output[[tag]], `[`, 1L))
    output[[tag]] <- as.data.frame(output[[tag]], row.names = basename(file),
                                   stringsAsFactors = FALSE)
  }

  # Coerce numeric_version information in Header
  cols <- intersect(c("FileVersion", "SoftwareVersion"),
                    colnames(output[["Header"]]))
  if (length(cols)) {
    output[["Header"]][, cols] <- lapply(output[["Header"]][, cols],
                                         numeric_version)
  }

  # Coerce date data type in Sample_Attributes
  if ("Date" %in% colnames(output[["Sample_Attributes"]])) {
    output[["Sample_Attributes"]][["Date"]] <-
      as.Date(output[["Sample_Attributes"]][["Date"]], format = "%Y%m%d")
  }

  # Coerce numeric data in Lane Attributes
  cols <- intersect(c("ID", "FovCount", "FovCounted", "StagePosition"),
                    colnames(output[["Lane_Attributes"]]))
  if (length(cols)) {
    output[["Lane_Attributes"]][, cols] <-
      lapply(output[["Lane_Attributes"]][, cols], as.integer)
  }
  if ("BindingDensity" %in% colnames(output[["Lane_Attributes"]])) {
    output[["Lane_Attributes"]][["BindingDensity"]] <-
      as.numeric(output[["Lane_Attributes"]][["BindingDensity"]])
  }

  # Rename ID columns
  for (ptag in c("Sample", "Lane")) {
    tag <- sprintf("%s_Attributes", ptag)
    names(output[[tag]])[names(output[[tag]]) == "ID"] <- sprintf("%sID", ptag)
  }

  # Convert Code_Summary to data.frame object
  if (output[["Code_Summary"]][1L] != "CodeClass,Name,Accession,Count")
    stop("Code_Summary section header is not \"CodeClass,Name,Accession,Count\"")
  output[["Code_Summary"]][1L] <- "CodeClass,GeneName,Accession,Count"
  output[["Code_Summary"]] <- paste(output[["Code_Summary"]], collapse = "\n")
  output[["Code_Summary"]] <-
    read.csv(textConnection(output[["Code_Summary"]]),
             colClasses = c(CodeClass = "character", GeneName = "character",
                            Accession = "character", Count = "integer"))
  rownames(output[["Code_Summary"]]) <-
    sprintf("%s_%s_%s",
            output[["Code_Summary"]][["CodeClass"]],
            output[["Code_Summary"]][["GeneName"]],
            output[["Code_Summary"]][["Accession"]])

  output
}
