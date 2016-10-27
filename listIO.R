writelist<-function (x, filename = "data", append = FALSE, closefile = TRUE,
    outfile)
{
    if (!append) {
        outfile <- file(filename, "w")
        cat(file = outfile, append = append)
    }
    for (i in 1:length(x)) {
        if (!is.null(x[[i]])) {
            switch(data.class(x[[i]]), matrix = write.table(x[[i]],
                sep = "\t", file = outfile, append = TRUE), table = if (!is.null(names(x[[i]]))) {
                write.table(rbind(names(x[[i]]), x[[i]]), file = outfile,
                  append = TRUE, row.names = FALSE, col.names = FALSE,
                  sep = "\t")
            }
            else {
                write(x[[i]], file = outfile, append = TRUE)
            }, list = write.list(x[[i]], outfile = outfile, append = TRUE,
                closefile = FALSE), if (!is.null(names(x[[i]]))) {
                write.table(rbind(names(x[[i]]), x[[i]]), file = outfile,
                  append = TRUE, row.names = FALSE, col.names = FALSE,
                  sep = "\t")
            }
            else {
                write(x[[i]], ncolumns = length(x[[i]]), file = outfile,
                  append = TRUE)
            })
        }
    }
    if (closefile)
        close(outfile)
}

readlist<-function(listfile) {

 t=strsplit(readLines(listfile), split="\t")
 l=lapply(t, function(x){as.numeric(unlist(strsplit(x, split=" ")))})
 return(l)

}