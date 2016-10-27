library("org.Mm.eg.db")
library("GO.db")
library("annotate")
library("genefilter")
library("GOstats")
library("RColorBrewer")
library("xtable")
library("Rgraphviz")
library("KEGG.db")

setwd("/Users/StephyTse/Documents/DrugAbuse/addictive_drugs_effect/processeddata/NaiveSalineHeroinMorphine/Difgenes")
geneID_all <- scan("all_gene_geneID.txt")
onlySalineID <- scan("/Users/StephyTse/Documents/DrugAbuse/addictive_drugs_effect/processeddata/NaiveSalineHeroinMorphine/Difgenes/naiveHeroin/naiveHeroin_onlyHeroin_geneID.txt")

# Get the entrez gene identifiers that are mapped to a GO ID
x <- org.Mm.egGO
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes]) 

#find our geneID in geneID_all
inds <- which(names(xx) %in% geneID_all) 
gene_all<-names(xx[inds])

#same process for onlySaline
inds <- which(gene_all %in% onlySalineID) 
onlySaline<-gene_all[inds]

#standard HyperGeo
hgCutoff <- 1
#0.01
#set hgCutoff <-1
#
params <- new("GOHyperGParams",
              geneIds=onlySaline,
              universeGeneIds=gene_all,
              annotation="org.Mm.eg.db",
              ontology="BP",
              pvalueCutoff=hgCutoff,
              conditional=FALSE,
              testDirection="over")
GO_hgOver <- hyperGTest(params)
GO_hgOver
summary(GO_hgOver)
#htmlReport(GO_hgOver, file="/Users/StephyTse/Documents/DrugAbuse/addictive_drugs_effect/processeddata/NaiveSalineHeroinMorphine/Difgenes/naiveHeroin/naive/GO_hgOver_p_0_01.html")

#conditional GOHyperGeo
paramsCond <- params
conditional(paramsCond) <- TRUE
GO_hgCondOver <- hyperGTest(paramsCond)
GO_hgCondOver
GO_hgCondOver.sum=summary(GO_hgCondOver)
#htmlReport(GO_hgCondOver, file="/Users/StephyTse/Documents/DrugAbuse/addictive_drugs_effect/processeddata/NaiveSalineHeroinMorphine/Difgenes/naiveHeroin/naive/GO_hgCondOver_p_0_01.html")



if (nrow(GO_hgCondOver.sum)>0) {
    
    hgCondOver.sigGO=GO_hgCondOver.sum$GOBPID
    
    hgCondOver.symbol.all=NULL
    for (hgCondOver.sigGO1 in hgCondOver.sigGO) {
        #hgCondOver.sigGO1= hgCondOver.sigGO[31]
        cat(hgCondOver.sigGO1, "\n")
        hgCondOver1.geneid=geneIdsByCategory(GO_hgCondOver, catids = hgCondOver.sigGO1)[[1]]
        
        ## you need to revise the code below -- find a R package which can map entrez gene ids to
        ## gene symbols
        #hgCondOver1.affyid=affyid[match(hgCondOver1.geneid, entrezIds)]
        hgCondOver1.symbol=unlist(mget(hgCondOver1.geneid, org.Mm.egSYMBOL))
        names(hgCondOver1.symbol)=NULL
        
        if (length(hgCondOver1.symbol)==0){
            #catch case where no geneids found
            hgCondOver.symbol.all=rbind(hgCondOver.symbol.all, "no geneids found")
        }else{
            hgCondOver.symbol.all=rbind(hgCondOver.symbol.all, paste(hgCondOver1.symbol, collapse=", "))
        }
        
    }
    
    hgCondOver.sum.all=cbind(GO_hgCondOver.sum, hgCondOver.symbol.all)
    names(hgCondOver.sum.all)=c(names(GO_hgCondOver.sum), "Genes")
    summary.file=file.path("/Users/StephyTse/Documents/DrugAbuse/addictive_drugs_effect/processeddata/NaiveSalineHeroinMorphine/Difgenes/naiveHeroin/heroin", "hg_CondOver_sum_all_cutoff1.txt")
    write.table(hgCondOver.sum.all, file=summary.file, row.names=F, col.names=T, sep="\t", quote=F);
    summarygo.file=file.path("/Users/StephyTse/Documents/DrugAbuse/addictive_drugs_effect/processeddata/NaiveSalineHeroinMorphine/Difgenes/naiveHeroin/heroin", "hgCondOver_sum_term_cutoff1.txt")
    write.table(GO_hgCondOver.sum$Term, file=summarygo.file, row.names=F, col.names=F, sep="\t", quote=F);
    
    rm(GO_hgCondOver.sum, hgCondOver.sum.all)
    
}





#plotFuns
#step1_getMaxConnCompGraph_step2_coloredGoPlot <- function(hgOver, hgCondOver) {
	#step1

    #uGoDagRev <- ugraph(goDag(hgOver))
#   sigNodes <- sigCategories(hgOver)
    #adjNodes <- unlist(adj(uGoDagRev, sigNodes))
#   adjNodes <- unlist(adj(goDag(hgOver), sigNodes))
# displayNodes <- unique(c(sigNodes, adjNodes))
#   displayGraph <- subGraph(displayNodes, goDag(hgOver))
#   cc <- connComp(displayGraph)
#   ccSizes <- listLen(cc)
#   ccMaxIdx <- which(ccSizes == max(ccSizes))
#   ccMaxGraph <- subGraph(cc[[ccMaxIdx]], displayGraph)
#   ccMaxGraph

    #step2

#   nodeColors <- sapply(nodes(ccMaxGraph),
#                        function(n) {
#                            if (n %in% sigCategories(hgCondOver))
#                              "dark red"
#                            else if (n %in% sigCategories(hgOver))
#                              "pink"
#                            else
#                              "gray"
#                        })
#   nattr <- makeNodeAttrs(ccMaxGraph,
#                          label=nodes(ccMaxGraph),
#                          shape="ellipse",
#                          fillcolor=nodeColors,
#                          fixedsize=FALSE)
#   plot(ccMaxGraph, nodeAttrs=nattr)
#}
#step1_getMaxConnCompGraph_step2_coloredGoPlot(GO_hgOver, hgCondOver)

#KEGGHyperGeo
#hgCutoff <- 0.05
#params <- new("KEGGHyperGParams",
#             geneIds=onlySaline,
#             universeGeneIds=gene_all,
#             annotation="org.Mm.eg.db",
#             pvalueCutoff=hgCutoff,
#             testDirection="over")
#KEGG_hgOver <- hyperGTest(params)
#KEGG_hgOver
#summary(KEGG_hgOver)
#htmlReport(KEGG_hgOver, file="/Users/StephyTse/Documents/DrugAbuse/addictive_drugs_effect/processeddata/NSH/Difgenes/naivesaline/onlysaline/KEGG_hgOver.html")
