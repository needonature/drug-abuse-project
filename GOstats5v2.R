# nohup R --vanilla < GOstats3.R & ##run in back
# source("GOstats5.R")

#rm(list=ls())

library("KEGG.db")
library("Category")
#library("org.Hs.eg.db")########## but we use mouse 
library("org.Mm.eg.db")
library("GO.db")
library("annotate")
library("genefilter")
library("GOstats")
library("RColorBrewer")
library("xtable")

setwd("C:/Users/thinkpad/Desktop/drug abuse project")
# input #
resdir = "GO_output"
inputDir = "intput/"
#filePaths <- 
list.files(path=inputDir)#, pattern='*.txt')

# Read in  entrezIds for common universe file #
entrezUniverse <- read.delim("all_gene_geneID.txt", stringsAsFactors = F)
entrezUniverseIds <- unique(as.character(entrezUniverse[,1]))
# Remove universe set with no GO annotation #
GOmasterSet <- org.Mm.egGO
# Get the entrez gene identifiers that are mapped to a GO ID
mapped_genes <- mappedkeys(GOmasterSet)
# Convert to a list
GOmappedSet <- as.list(GOmasterSet[mapped_genes])
geneUniverse = GOmappedSet[entrezUniverseIds]
#Remove null entries
geneUniverseF =  names(geneUniverse[!geneUniverse=="NULL"])

###for (fileIndx in seq(1,length(filePaths))){

	# Read in entrezIDs from file
	geneCorFile <- filePaths[fileIndx]
	cat('Reading from ', fileIndx,'/',length(filePaths),': ', geneCorFile,'\n');
	rawDF <- read.delim( paste0(inputDir, geneCorFile), stringsAsFactors = F)
	entrezIdsRaw <- rawDF$EntrezGeneID[which(!is.na(rawDF$EntrezGeneID))]
	entrezIds <- unique(as.character(entrezIdsRaw)) 
	
	
	# remove probe sets for which we have no GO annotation #
	uniqGenes = GOmappedSet[entrezIds]
	#Remove null entries
	uniqGenesF = names(uniqGenes[!uniqGenes=="NULL"])

  
  directions=c("over", "under")
  hgCutoff <- 0.005

	
###  for (direction in directions) {
    direction=directions[1]

      cat("Calculating", direction, "-represented genes", "\n")
      params <- new("GOHyperGParams", geneIds = uniqGenesF,
          universeGeneIds = geneUniverseF, annotation = "org.Mm.eg.db",
          ontology = "BP", pvalueCutoff = hgCutoff, conditional = FALSE,
          testDirection = direction)

      paramsCond <- params
      conditional(paramsCond) <- TRUE

      #Outputs and Result Summarization#
      tryCatch({
          hgCondOver <- hyperGTest(paramsCond)

          # write to files #

          if (length(hgCondOver@geneIds)>0) {

          #resdir=file.path(datadir1.fpath, "results")
          if (!file.exists(resdir)) {
          dir.create(resdir)
      }
      celltype=gsub("BAL_sigcorBALgenes_", "",geneCorFile)
      celltype=gsub(".txt", "", celltype)
      
	    html.file=file.path(resdir, paste("GOCond", direction, '_',celltype, ".html", sep=""))
      htmlReport(hgCondOver, file = html.file, summary.args=list("htmlLinks"=TRUE))

      hgCondOver.sum=summary(hgCondOver)

## Siwei and Xin, you need the code below to get genes associated with enriched GO terms ##


    if (nrow(hgCondOver.sum)>0) {

      hgCondOver.sigGO=hgCondOver.sum$GOBPID

      hgCondOver.symbol.all=NULL
      for (hgCondOver.sigGO1 in hgCondOver.sigGO) {
      # hgCondOver.sigGO1= hgCondOver.sigGO[31]
      cat(hgCondOver.sigGO1, "\n")
  	  hgCondOver1.geneid=geneIdsByCategory(hgCondOver, catids = hgCondOver.sigGO1)[[1]]

    	## you need to revise the code below -- find a R package which can map entrez gene ids to gene symbols 
    	#hgCondOver1.affyid=affyid[match(hgCondOver1.geneid, entrezIds)]
      hgCondOver1.symbol=unlist(mget(hgCondOver1.geneid, org.Hs.egSYMBOL))
      names(hgCondOver1.symbol)=NULL
    	
    	if (length(hgCondOver1.symbol)==0){
        	  #catch case where no geneids found 
        		hgCondOver.symbol.all=rbind(hgCondOver.symbol.all, "no geneids found")
      }else{
        	 hgCondOver.symbol.all=rbind(hgCondOver.symbol.all, paste(hgCondOver1.symbol, sep=", "))
    	}
	
      }

      hgCondOver.sum.all=cbind(hgCondOver.sum, hgCondOver.symbol.all)
      names(hgCondOver.sum.all)=c(names(hgCondOver.sum), "Genes")
      summary.file=file.path(resdir, paste("GoCond", direction,"sum_", celltype, ".txt", sep=""))
      write.table(hgCondOver.sum.all, file=summary.file, row.names=F, col.names=T, sep="\t", quote=F);
      summarygo.file=file.path(resdir, paste("GoCond", direction,"sumBP_", celltype, ".txt", sep=""))
      write.table(hgCondOver.sum$Term, file=summarygo.file, row.names=F, col.names=F, sep="\t", quote=F);

      rm(hgCondOver.sum, hgCondOver.sum.all)

    } # if (nrow(hgCondOver.sum)>0)

    rm(hgCondOver)

   } # end if (length ...

},error=function(err){
	cat(paste("error:  ",err," Will not write file.\n"))
})
###} # end direction
###}
 