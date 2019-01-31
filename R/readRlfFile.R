readRlfFile <-
function(file)
{
  # Read data from Reporter Library File (RLF) file
  lines <- trimws(readLines(file))

  # Split data by tags
  locs <- match(c("[Header]", "[Content]"), lines)
  contentDim <- lines[(locs[2L] + 2L):(locs[2L] + 1L)]
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
                 SpotsPerBarcode = as.integer(sections[["Header"]][["NSpot"]]),
                 BasePairsPerSpot = as.integer(sections[["Header"]][["NBasePair"]]),
                 BackboneType = sections[["Header"]][["Backbone"]],
                 CodeClassCount = as.integer(sections[["Header"]][["ClassCount"]]))

  # Create Classes
  k <- header[["CodeClassCount"]]
  getTagValues <- function(tag) sections[["Header"]][sprintf(tag, 0L:(k-1L))]
  classes <-
    data.frame(Classification = as.integer(getTagValues("ClassKey%d")),
               CodeClass = getTagValues("ClassName%d"),
               CodeClassActive = as.integer(getTagValues("ClassActive%d")),
               CodeClassDate = as.Date(getTagValues("ClassDate%d"),
                                       format = "%Y%m%d"),
               CodeClassSource = getTagValues("ClassSource%d"),
               CodeClassPreparer = getTagValues("ClassPreparer%d"),
               stringsAsFactors = FALSE)

  # Parse Content
  nrecords <- length(sections[["Content"]]) - 1L
  if (contentDim[1L] != sprintf("RecordCount=%d", nrecords))
    stop(sprintf("Content section RecordCount must be %d", nrecords))
  if (contentDim[2L] != "ColumnCount=8")
    stop("Content section ColumnCount must be 8")
  if (sections[["Content"]][[1L]][1L] != "Classification,TargetSeq,BarCode,GeneName,ProbeID,Species,Accession,Comments")
    stop("Content section header is not \"Classification,TargetSeq,BarCode,GeneName,ProbeID,Species,Accession,Comments\"")
  sections[["Content"]][[1L]][1L] <- "Classification,TargetSeq,Barcode,GeneName,ProbeID,Species,Accession,BarcodeComments"
  rn <- names(sections[["Content"]])[-1L]
  output <-
    read.csv(textConnection(paste(sections[["Content"]], collapse = "\n")),
             row.names = rn,
             colClasses = c(Classification = "integer", TargetSeq = "character",
                            Barcode = "character", GeneName = "character",
                            ProbeID = "character", Species = "character",
                            Accession = "character",
                            BarcodeComments = "character"))

  # Merge classes information with output
  output[["RowNames"]] <- rn
  output <- merge(output, classes, by = "Classification", all.x = TRUE,
                   sort = FALSE)
  rownames(output) <- output[["RowNames"]]
  output[["Classification"]] <- output[["RowNames"]] <- NULL
  output <- output[rn, c("CodeClass", "GeneName", "Accession", "TargetSeq",
                         "Barcode", "ProbeID", "Species",  "BarcodeComments",
                         "CodeClassActive", "CodeClassDate", "CodeClassSource",
                         "CodeClassPreparer")]

  # Coerce output to DataFrame
  output <- DataFrame(output)

  # Convert sequences to XStringSet columns
  output[["TargetSeq"]] <-
    DNAStringSet(ifelse(is.na(output[["TargetSeq"]]), ".",
                        output[["TargetSeq"]]))
  output[["Barcode"]] <-
    BStringSet(ifelse(is.na(output[["Barcode"]]), ".", output[["Barcode"]]))

  # Move constant columns to header
  if (all(output[["CodeClassDate"]] == header[["RlfFileDate"]], na.rm = TRUE)) {
    output[["CodeClassDate"]] <- NULL
  }
  for (j in c("CodeClassSource", "CodeClassPreparer")) {
    if (length(value <- unique(output[[j]])) == 1L) {
      header[[j]] <- value
      output[[j]] <- NULL
    }
  }

  # Add header to metadata
  metadata(output) <- header

  output
}
