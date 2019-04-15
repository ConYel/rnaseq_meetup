
UpdatePrefix <- function(prefix, ...)
{
    dots <- unlist(list(...))
    if(length(dots) != 0)
    {
        for (str in dots)
        {
            # str <- gsub(pattern=".", replacement="_", str)
            prefix <- paste(prefix, str, sep=" ")
        }

    } else {
        stop("provide a string to append to ", new.prefix)
    }
    return(prefix)
}

UpdateFolderPath <- function(path, ...)
{
    dots <- unlist(list(...))
    if(length(dots) != 0)
    {
        for (str in dots)
        {
            str <- gsub(pattern=" ", replacement="_", str)
            path <- file.path(path, str)
        }
    } else {
        stop("provide a string to append to ", path)
    }
    dir.create(path, recursive=TRUE, showWarnings=FALSE)
    message("Recursively created ", path, " on disk")
    return(path)
}

UpdateFilename <- function(filename, ..., extension=NULL)
{
    dots <- unlist(list(...))
    print(dots)
    filename <- gsub(pattern=" ", replacement="_", x=filename)
    if(length(dots) != 0)
    {
        for (str in dots) filename <- paste(filename, str, sep="_")
    } else {
        stop("provide a string to append to ", filename)
    }
    if(!is.null(extension))
    {
        filename <- paste0(filename, ".", extension)
    }
    return(filename)
}

Venn2de <-
  function(x, y, label1, label2, title="VENN", plot.dir="./", conversion.map=NULL)
{

    require(limma)

    out.path <- UpdateFolderPath(plot.dir, "venn2")

    a15 = x

    b15 = y

    c15 <- intersect(a15, b15)     #common gene names

    ab <- setdiff(a15, b15)

    ba <- setdiff(b15, a15)


    Lists <- list(a15, b15)  #put the word vectors into a list to supply lapply
    Lists <- lapply(Lists, function(x) as.character(unlist(x)))
    items <- sort(unique(unlist(Lists)))   #put in alphabetical order
    MAT <- matrix(rep(0, length(items)*length(Lists)), ncol=2)  #make a matrix of 0's
    names <- c(label1,label2)
    colnames(MAT) <- names
    rownames(MAT) <- items
    lapply(seq_along(Lists), function(i) {   #fill the matrix
      MAT[items %in% Lists[[i]], i] <<- table(Lists[[i]])
    })

    outputName=paste(label1,"_",label2,"_genes_in_intersection.txt",sep="")

    outputName2=paste("genes_in_",label1,"_not_in_",label2,".txt",sep="")

    outputName3=paste("genes_in_",label2,"_not_in_",label1,".txt",sep="")

    outputName <- file.path(out.path, outputName)
    outputName2 <- file.path(out.path, outputName2)
    outputName3 <- file.path(out.path, outputName3)
    if(!is.null(conversion.map)){
      c15 <- CreateConvertedDataframe(c15, conversion.map = conversion.map)
      ab <- CreateConvertedDataframe(ab, conversion.map = conversion.map)
      ba <- CreateConvertedDataframe(ba, conversion.map = conversion.map)
    }

    write.table(c15, file = outputName , quote=FALSE, sep="\t", row.names=FALSE)

    write.table(ab, file = outputName2 , quote=FALSE, sep="\t", row.names=FALSE)

    write.table(ba, file = outputName3 , quote=FALSE, sep="\t", row.names=FALSE)


    outputName=file.path(out.path, paste(label1,"_",label2,"_VennDiagramDE.pdf",sep=""))

    # dev.new()
    limma::vennDiagram(MAT, circle.col= c("red","green"), main=title)

    # dev.print(device = pdf, file=outputName, width=10, height=10)
    # dev.off()
    # print(c15)
    return(c15)

  }
