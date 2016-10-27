library("org.Mm.eg.db")
library("GO.db")
library("annotate")
library("genefilter")
library("GOstats")
library("RColorBrewer")
library("xtable")
library("Rgraphviz")
library("base")
#library("KEGG.db")
#biocLite("")

#setwd("C:\\Users\\thinkpad\\Desktop\\drug abuse project")
setwd("/Users/guest1/xin")
geneID_all <- scan("all_gene_geneID.txt")
naiveCocaine_onlycocaine_geneID <- scan("naiveCocaine_onlycocaine_geneID.txt")

# Get the entrez gene identifiers that are mapped to a GO ID
x <- org.Mm.egGO
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes]) 

#find our geneID in geneID_all
inds <- which(names(xx) %in% geneID_all) 
gene_all<-names(xx[inds])

#same process for onlySaline
inds <- which(gene_all %in% naiveCocaine_onlycocaine_geneID) 
onlySaline<-gene_all[inds]

#standard HyperGeo
hgCutoff <- 0.01  #when we set 1, we get all the GO term's p-value
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
#htmlReport(GO_hgOver, file="/Users/StephyTse/Documents/DrugAbuse/addictive_drugs_effect/processeddata/NSH/Difgenes/salineHeroin/onlysaline/GO_hgOver_p_0_01.html")

#conditional GOHyperGeo
paramsCond <- params
conditional(paramsCond) <- TRUE
GO_hgCondOver <- hyperGTest(paramsCond)
GO_hgCondOver
GO_hgCondOver.sum=summary(GO_hgCondOver)
#htmlReport(GO_hgCondOver, file="/Users/guest1/xin/GO_hgCondOver__naiveCocaine_onlycocaine.html")

#get genes associated with enriched GO terms
if (nrow(GO_hgCondOver.sum)>0) {    
    hgCondOver.sigGO=GO_hgCondOver.sum$GOBPID
    hgCondOver.symbol.all=NULL
    #hgCondOver.sigGO1=hgCondOver.sigGO[1]
    for (hgCondOver.sigGO1 in hgCondOver.sigGO) {
        #hgCondOver.sigGO1= hgCondOver.sigGO[1]
        cat(hgCondOver.sigGO1, "\n")
        hgCondOver1.geneid=geneIdsByCategory(GO_hgCondOver, catids = hgCondOver.sigGO1)[[1]]
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
    summary.file=file.path("/Users/guest1/xin/", "hg_CondOver_sum_all_GO_naiveCocaine_onlycocaine.txt")
    write.table(hgCondOver.sum.all, file=summary.file, row.names=F, col.names=T, sep="\t", quote=F);
    summarygo.file=file.path("/Users/guest1/xin/", "hgCondOver_sum_term_GO_naiveCocaine_onlycocaine.txt")
    write.table(GO_hgCondOver.sum$Term, file=summarygo.file, row.names=F, col.names=F, sep="\t", quote=F);
    
    #rm(GO_hgCondOver.sum, hgCondOver.sum.all)
}

#for adjusting p values, first find all go term
pvalueCutoff(paramsCond) <- 1
GO_hgCondOver_allGoTerm <- hyperGTest(paramsCond)
GO_hgCondOver.sum_adpv<-summary(GO_hgCondOver_allGoTerm)
summary.file=file.path("/Users/guest1/xin/", "hg_CondOver_sum_all_allGOterm_GO_naiveCocaine_onlycocaine.txt")
write.table(GO_hgCondOver.sum_adpv, file=summary.file, row.names=F, col.names=T, sep="\t", quote=F);
summarygo.file=file.path("/Users/guest1/xin/", "hgCondOver_sum_term_allGOterm_GO_naiveCocaine_onlycocaine.txt")
write.table(GO_hgCondOver.sum_adpv$Term, file=summarygo.file, row.names=F, col.names=F, sep="\t", quote=F);
rawp=GO_hgCondOver.sum_adpv$Pvalue

#adjusting p values
adjp=p.adjust(rawp, method="BH")
GO_hgCondOver.sum_adpv$Pvalue<-adjp

#get genes associated with enriched GO terms
if (nrow(GO_hgCondOver.sum_adpv)>0) {    
    hgCondOver.sigGO=GO_hgCondOver.sum_adpv$GOBPID
    hgCondOver.symbol.all=NULL
    #hgCondOver.sigGO1=hgCondOver.sigGO[1]
    for (hgCondOver.sigGO1 in hgCondOver.sigGO) {
        #hgCondOver.sigGO1= hgCondOver.sigGO[1]
        cat(hgCondOver.sigGO1, "\n")
        hgCondOver1.geneid=geneIdsByCategory(GO_hgCondOver_allGoTerm, catids = hgCondOver.sigGO1)[[1]]
        hgCondOver1.symbol=unlist(mget(hgCondOver1.geneid, org.Mm.egSYMBOL))
        names(hgCondOver1.symbol)=NULL

        if (length(hgCondOver1.symbol)==0){
            #catch case where no geneids found
            hgCondOver.symbol.all=rbind(hgCondOver.symbol.all, "no geneids found")
        }else{
            hgCondOver.symbol.all=rbind(hgCondOver.symbol.all, paste(hgCondOver1.symbol, collapse=", "))
        }       
    }
    
    hgCondOver.sum.all=cbind(GO_hgCondOver.sum_adpv, hgCondOver.symbol.all)
    names(hgCondOver.sum.all)=c(names(GO_hgCondOver.sum_adpv), "Genes")
    summary.file=file.path("/Users/guest1/xin/", "hg_CondOver_sum_all_adpv_GO_naiveCocaine_onlycocaine.txt")
    write.table(hgCondOver.sum.all, file=summary.file, row.names=F, col.names=T, sep="\t", quote=F);
    summarygo.file=file.path("/Users/guest1/xin/", "hgCondOver_sum_term_adpv_GO_naiveCocaine_onlycocaine.txt")
    write.table(GO_hgCondOver.sum_adpv$Term, file=summarygo.file, row.names=F, col.names=F, sep="\t", quote=F);

    #rm(GO_hgCondOver.sum_adpv, hgCondOver.sum.all)
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
