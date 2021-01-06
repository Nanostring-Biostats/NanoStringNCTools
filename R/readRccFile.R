readRccFile <- function(file) {
    lines <- trimws(readLines(file))
    tags <- names(.rccMetadata[["schema"]])
    output <- lapply(tags, function(tag) {
        bounds <- charmatch(sprintf(c("<%s>", "</%s>"), tag), lines)
        if (anyNA(bounds) || bounds[1L] + 1L >= bounds[2L]) 
            lines[integer(0)]
        else lines[(bounds[1L] + 1L):(bounds[2L] - 1L)]
    })
    names(output) <- tags
    for (tag in c("Header", "Sample_Attributes", "Lane_Attributes")) {
        while (length(bad <- grep(",", output[[tag]], invert = TRUE)) > 0L) {
            bad <- bad[1L]
            if (bad == 1L) 
                stop(sprintf("%s section has malformed first line", tag))
            fixed <- output[[tag]]
            fixed[bad - 1L] <- sprintf("%s %s", fixed[bad - 1L], fixed[bad])
            output[[tag]] <- fixed[-bad]
        }
        output[[tag]] <- strsplit(output[[tag]], split = ",")
        output[[tag]] <- structure(lapply(output[[tag]], function(x) if (length(x) == 1L) 
            ""
        else x[2L]), names = lapply(output[[tag]], `[`, 1L), class = "data.frame", row.names = basename(file))
    }
    metadata <- .rccMetadata[["schema"]][["Header"]]
    cols <- c("FileVersion", "SoftwareVersion")
    if (!(all(cols %in% colnames(output[["Header"]])))) 
        stop("Header section must contain \"FileVersion\" and \"SoftwareVersion\"")
    output[["Header"]][, cols] <- lapply(output[["Header"]][, cols], numeric_version)
    fileVersion <- output[["Header"]][1L, "FileVersion"]
    if (!(fileVersion %in% numeric_version(c("1.7", "2.0")))) 
        stop("\"FileVersion\" in Header section must be either 1.7 or 2.0")
    for (section in c("Header", "Sample_Attributes", "Lane_Attributes")) {
        valid <- .validRccSchema(output[[section]], fileVersion, section)
        if (!isTRUE(valid)) 
            stop(valid)
    }
    output[["Sample_Attributes"]][["Date"]] <- as.Date(output[["Sample_Attributes"]][["Date"]], 
        format = "%Y%m%d")
    cols <- c("ID", "FovCount", "FovCounted", "StagePosition")
    output[["Lane_Attributes"]][, cols] <- lapply(output[["Lane_Attributes"]][, cols], 
        as.integer)
    output[["Lane_Attributes"]][["BindingDensity"]] <- as.numeric(output[["Lane_Attributes"]][["BindingDensity"]])
    for (ptag in c("Sample", "Lane")) {
        tag <- sprintf("%s_Attributes", ptag)
        names(output[[tag]])[names(output[[tag]]) == "ID"] <- sprintf("%sID", ptag)
    }
    for (col in c("Owner", "Comments", "Date")) {
        names(output[["Sample_Attributes"]])[names(output[["Sample_Attributes"]]) == col] <- sprintf("Sample%s", 
            col)
    }
    if (fileVersion == numeric_version("1.7")) {
        output[["Header"]][["SystemType"]] <- "Gen2"
        output[["Sample_Attributes"]][["AssayType"]] <- NA_character_
    }
    if (output[["Code_Summary"]][1L] != "CodeClass,Name,Accession,Count") 
        stop("Code_Summary section header is not \"CodeClass,Name,Accession,Count\"")
    output[["Code_Summary"]][1L] <- "CodeClass,GeneName,Accession,Count"
    output[["Code_Summary"]] <- paste(output[["Code_Summary"]], collapse = "\n")
    output[["Code_Summary"]] <- read.csv(textConnection(output[["Code_Summary"]]), colClasses = c(CodeClass = "character", 
        GeneName = "character", Accession = "character", Count = "numeric"))
    output[["Code_Summary"]][["Count"]] <- as.integer(round(output[["Code_Summary"]][["Count"]]))
    rn <- sprintf("%s_%s_%s", output[["Code_Summary"]][["CodeClass"]], output[["Code_Summary"]][["GeneName"]], 
        output[["Code_Summary"]][["Accession"]])
    if ((ndups <- anyDuplicated(rn)) > 0L) {
        warning(sprintf("removed %d rows from \"Code_Summary\" due to duplicate rownames", 
            ndups))
        ok <- which(!duplicated(rn, fromLast = FALSE) & !duplicated(rn, fromLast = TRUE))
        rn <- rn[ok]
        output[["Code_Summary"]] <- output[["Code_Summary"]][ok, , drop = FALSE]
    }
    rownames(output[["Code_Summary"]]) <- rn
    output
}
