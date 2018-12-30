bgThresh <- function(x, STATS) sweep(x, 2L, STATS, FUN = "pmax")
bgSubtr  <- function(x, STATS) sweep(x, 2L, STATS, FUN = "-")

fgScale <- function(x, STATS) sweep(x, 2L, STATS / mean(STATS), FUN = "/")
