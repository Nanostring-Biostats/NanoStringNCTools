ggplot.NanoStringRccSet <- function(data, mapping = aes(), ..., extradata = NULL, elt = "exprs", 
    tooltip_digits = 4L, environment = parent.frame()) {
    if (length(mapping) == 0L) {
        mapping <- design(data)
        if (is.null(mapping)) 
            stop("\"mapping\" argument is missing")
    }
    df <- munge(data, mapping = mapping, extradata = extradata, elt = elt)
    if ("tooltip" %in% names(mapping)) {
        tooltip <- as.character(mapping[["tooltip"]][[2L]])
        for (j in c("x", "y")) {
            if (j %in% names(mapping)) {
                mf <- model.frame(mapping[[j]], df)
                df[[tooltip]] <- sprintf("%s | %s&nbsp;=&nbsp;%s", df[[tooltip]], names(mf)[1L], 
                  signif(mf[[1L]], digits = tooltip_digits))
            }
        }
    }
    ggplot(df, mapping, ..., environment = environment)
}
geom_beeswarm_interactive <- function(mapping = NULL, data = NULL, priority = c("ascending", 
    "descending", "density", "random", "none"), cex = 1, groupOnX = NULL, dodge.width = 0, 
    stat = "identity", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, ...) {
    position <- position_beeswarm(priority = priority, cex = cex, groupOnX = groupOnX, 
        dodge.width = dodge.width)
    layer(data = data, mapping = mapping, stat = stat, geom = GeomInteractivePoint, position = position, 
        show.legend = show.legend, inherit.aes = inherit.aes, params = list(na.rm = na.rm, 
            ...))
}
update_geom_params <- function(geom, new, old = aes()) {
    if (geom %in% names(new)) {
        if ("color" %in% names(new[[geom]])) {
            new[[geom]][["colour"]] <- new[[geom]][["color"]]
            new[[geom]][["color"]] <- NULL
        }
    }
    else {
        new[[geom]] <- aes()
    }
    if (length(old) > 0L) {
        args <- setdiff(names(old), c("tooltip", "onclick", "data_id"))
        new[[geom]] <- new[[geom]][names(new[[geom]]) %in% args]
        new[[geom]] <- c(new[[geom]], old[setdiff(args, names(new[[geom]]))])
        oldClass(new[[geom]]) <- "uneval"
    }
    new
}
format_percent <- function(x) {
    sprintf("%s%%", format(100 * x, digits = 2L))
}
