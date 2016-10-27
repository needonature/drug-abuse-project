library("org.Mm.eg.db")
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
x <- org.Mm.egPATH 
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes]) 

#find our geneID in geneID_all
inds <- which(names(xx) %in% geneID_all) 
gene_all<-names(xx[inds])

#same process for onlySaline
inds <- which(gene_all %in% onlySalineID) 
onlySaline<-gene_all[inds]

hgCutoff <- 1
params <- new("KEGGHyperGParams",
             geneIds=onlySaline,
             universeGeneIds=gene_all,
             annotation="org.Mm.eg.db",
             pvalueCutoff=hgCutoff,
             testDirection="over")
KEGG_hgOver <- hyperGTest(params)
KEGG_hgOver
Kegg_hgOver.sum=summary(KEGG_hgOver)

if (nrow(Kegg_hgOver.sum)>0){


    hgOver.sigkegg=Kegg_hgOver.sum$KEGGID
    
    hgOver.symbol.all=NULL
    for (hgOver.sigkegg1 in hgOver.sigkegg) {
        #hgCondOver.sigGO1= hgCondOver.sigGO[31]
        cat(hgOver.sigkegg1, "\n")
        hgOver1.geneid=geneIdsByCategory(KEGG_hgOver, catids = hgOver.sigkegg1)[[1]]
        
        ## you need to revise the code below -- find a R package which can map entrez gene ids to
        ## gene symbols
        #hgCondOver1.affyid=affyid[match(hgCondOver1.geneid, entrezIds)]
        hgOver1.symbol=unlist(mget(hgOver1.geneid, org.Mm.egSYMBOL))
        names(hgOver1.symbol)=NULL
        
        if (length(hgOver1.symbol)==0){
            #catch case where no geneids found
            hgOver.symbol.all=rbind(hgOver.symbol.all, "no geneids found")
        }else{
            hgOver.symbol.all=rbind(hgOver.symbol.all, paste(hgOver1.symbol, collapse=", "))
        }
        
    }
    
    hgOver.sum.all=cbind(Kegg_hgOver.sum, hgOver.symbol.all)
    names(hgOver.sum.all)=c(names(Kegg_hgOver.sum), "Genes")
    summary.file=file.path("/Users/StephyTse/Documents/DrugAbuse/addictive_drugs_effect/processeddata/NaiveSalineHeroinMorphine/Difgenes/naiveHeroin/heroin", "kegg_over_sum_all_cutoff1.txt")
    write.table(hgCondOver.sum.all, file=summary.file, row.names=F, col.names=T, sep="\t", quote=F);
    summarygo.file=file.path("/Users/StephyTse/Documents/DrugAbuse/addictive_drugs_effect/processeddata/NaiveSalineHeroinMorphine/Difgenes/naiveHeroin/heroin", "kegg_over_sum_term_cutoff1.txt")
    write.table(Kegg_hgOver.sum$Term, file=summarygo.file, row.names=F, col.names=F, sep="\t", quote=F);
    
    rm(Kegg_hgOver.sum, hgOver.sum.all)
    

}