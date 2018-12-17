readRlfFile <-
function(file)
{
  # Read data from Reporter Library File (RLF) file
  lines <- trimws(readLines(file))

  # Split data by tags
  locs <- match(c("[Header]", "[Content]"), lines)
  sections <-
    list("Header"  = lines[(locs[1L] + 1L):(locs[2L] - 1L)],
         "Content" = lines[(locs[2L] + 3L):length(lines)])
  sections <- lapply(sections, function(x)
    structure(sub(".*=(.*)", "\\1", x),
              names = sub("(.*)=.*", "\\1", x)))

  # Extract header
  header <- list(RlfFileVersion =
                   numeric_version(sections[["Header"]][["Version"]]),
                 RlfFileDate = as.Date(sections[["Header"]][["Date"]],
                                       format = "%Y%m%d"),
                 SpotsPerCode = as.integer(sections[["Header"]][["NSpot"]]),
                 BasePairsPerSpot = as.integer(sections[["Header"]][["NBasePair"]]),
                 BackboneType = sections[["Header"]][["Backbone"]],
                 CodeClassCount = as.integer(sections[["Header"]][["ClassCount"]]))

  # Create Classes
  k <- header[["CodeClassCount"]]
  getTagValues <- function(tag) sections[["Header"]][sprintf(tag, 0L:(k-1L))]
  classes <-
    data.frame(Classification = as.integer(getTagValues("ClassKey%d")),
               CodeClass = getTagValues("ClassName%d"),
               Active = as.integer(getTagValues("ClassActive%d")),
               Date = as.Date(getTagValues("ClassDate%d"), format = "%Y%m%d"),
               Source = getTagValues("ClassSource%d"),
               Preparer = getTagValues("ClassPreparer%d"),
               stringsAsFactors = FALSE)

  # Parse Content
  if (sections[["Content"]][[1L]][1L] != "Classification,TargetSeq,BarCode,GeneName,ProbeID,Species,Accession,Comments")
    stop("Content section header is not \"Classification,TargetSeq,BarCode,GeneName,ProbeID,Species,Accession,Comments\"")
  rn <- names(sections[["Content"]])[-1L]
  output <-
    read.csv(textConnection(paste(sections[["Content"]], collapse = "\n")),
             row.names = rn,
             colClasses = c(Classification = "integer", TargetSeq = "character",
                            BarCode = "character", GeneName = "character",
                            ProbeID = "character", Species = "character",
                            Accession = "character", Comments = "character"))

  # Merge classes information with output
  output[["RowNames"]] <- rn
  output <- merge(output, classes, by = "Classification", all.x = TRUE,
                   sort = FALSE)
  rownames(output) <- output[["RowNames"]]
  output[["Classification"]] <- output[["RowNames"]] <- NULL
  output <- output[rn, c("CodeClass", "GeneName", "Accession", "TargetSeq",
                         "BarCode", "ProbeID", "Species",  "Comments" ,
                         "Active", "Date", "Source", "Preparer")]

  # Coerce output to DataFrame
  output <- DataFrame(output)

  # Convert sequences to XStringSet columns
  if ("TargetSeq" %in% colnames(output)) {
    output[["TargetSeq"]] <-
      DNAStringSet(ifelse(is.na(output[["TargetSeq"]]), ".",
                          output[["TargetSeq"]]))
  }
  if ("BarCode" %in% colnames(output)) {
    output[["BarCode"]] <-
      BStringSet(ifelse(is.na(output[["BarCode"]]), ".", output[["BarCode"]]))
  }

  # Move constant columns to header
  if (all(output[["Date"]] == header[["Date"]], na.rm = TRUE)) {
    output[["Date"]] <- NULL
  }
  for (j in c("Source", "Preparer")) {
    if (length(value <- unique(output[[j]])) == 1L) {
      header[[j]] <- value
      output[[j]] <- NULL
    }
  }

  # Add header to metadata
  metadata(output) <- header

  output
}
