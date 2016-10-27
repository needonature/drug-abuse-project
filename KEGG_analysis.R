library("org.Mm.eg.db")
library("GO.db")
library("annotate")
library("genefilter")
library("GOstats")
library("RColorBrewer")
library("xtable")
library("Rgraphviz")
library("KEGG.db")
library("base")
#biocLite("")

#setwd("C:\\Users\\thinkpad\\Desktop\\drug abuse project")
setwd("/Users/guest1/xin")
geneID_all <- scan("all_gene_geneID.txt")
naiveCocaine_onlycocaine_geneID <- scan("naiveCocaine_onlycocaine_geneID.txt")

# Get the entrez gene identifiers that are mapped to a KEGG ID
x <- org.Mm.egPATH
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes]) 

#find our geneID in geneID_all
inds <- which(names(xx) %in% geneID_all) 
gene_all<-names(xx[inds])

#same process for onlySaline
inds <- which(gene_all %in% naiveCocaine_onlycocaine_geneID) 
onlySaline<-gene_all[inds]


hgCutoff <- 0.05
params <- new("KEGGHyperGParams",
            geneIds=onlySaline,
            universeGeneIds=gene_all,
            annotation="org.Mm.eg.db",
            pvalueCutoff=hgCutoff,
            testDirection="over")
KEGG_hgOver <- hyperGTest(params)
KEGG_hgOver
KEGG_hgOver.sum<-summary(KEGG_hgOver)
#htmlReport(KEGG_hgOver, file="/Users/StephyTse/Documents/DrugAbuse/addictive_drugs_effect/processeddata/NSH/Difgenes/naivesaline/onlysaline/KEGG_hgOver.html")


#get genes associated with enriched GO terms
if (nrow(KEGG_hgOver.sum)>0) {    
    hgCondOver.sigGO=KEGG_hgOver.sum$GOBPID
    hgCondOver.symbol.all=NULL
    #hgCondOver.sigGO1=hgCondOver.sigGO[1]
    for (hgCondOver.sigGO1 in hgCondOver.sigGO) {
        #hgCondOver.sigGO1= hgCondOver.sigGO[5]
        cat(hgCondOver.sigGO1, "\n")
        hgCondOver1.geneid=geneIdsByCategory(KEGG_hgOver, catids = hgCondOver.sigGO1)[[1]]
        hgCondOver1.symbol=unlist(mget(hgCondOver1.geneid, org.Mm.egSYMBOL))
        names(hgCondOver1.symbol)=NULL

        if (length(hgCondOver1.symbol)==0){
            #catch case where no geneids found
            hgCondOver.symbol.all=rbind(hgCondOver.symbol.all, "no geneids found")
        }else{
            hgCondOver.symbol.all=rbind(hgCondOver.symbol.all, paste(hgCondOver1.symbol, collapse=", "))
        }       
    }
    
    hgCondOver.sum.all=cbind(KEGG_hgOver.sum, hgCondOver.symbol.all)
    names(hgCondOver.sum.all)=c(names(KEGG_hgOver.sum), "Genes")
    summary.file=file.path("/Users/guest1/xin/", "hg_CondOver_sum_all_KEGG_naiveCocaine_onlycocaine.txt")
    write.table(hgCondOver.sum.all, file=summary.file, row.names=F, col.names=T, sep="\t", quote=F);
    summarygo.file=file.path("/Users/guest1/xin/", "hgCondOver_sum_term_KEGG_naiveCocaine_onlycocaine.txt")
    write.table(KEGG_hgOver.sum$Term, file=summarygo.file, row.names=F, col.names=F, sep="\t", quote=F);
    #rm(KEGG_hgOver.sum, hgCondOver.sum.all)
}

#for adjusting p values, first find all go term
pvalueCutoff(paramsCond) <- 1
KEGG_hgOver_allGoTerm <- hyperGTest(paramsCond)
KEGG_hgOver.sum_adpv<-summary(KEGG_hgOver_allGoTerm)
summary.file=file.path("/Users/guest1/xin/", "hg_CondOver_sum_all_allGOterm_KEGG_naiveCocaine_onlycocaine.txt")
write.table(GO_hgCondOver.sum_adpv, file=summary.file, row.names=F, col.names=T, sep="\t", quote=F);
summarygo.file=file.path("/Users/guest1/xin/", "hgCondOver_sum_term_allGOterm_KEGG_naiveCocaine_onlycocaine.txt")
write.table(GO_hgCondOver.sum_adpv$Term, file=summarygo.file, row.names=F, col.names=F, sep="\t", quote=F);
rawp=KEGG_hgOver.sum_adpv$Pvalue

#adjusting p values
adjp=p.adjust(rawp, method="hochberg")
KEGG_hgOver.sum_adpv$Pvalue<-adjp

#get genes associated with enriched GO terms
if (nrow(KEGG_hgOver.sum_adpv)>0) {    
    hgCondOver.sigGO=KEGG_hgOver.sum_adpv$GOBPID
    hgCondOver.symbol.all=NULL
    #hgCondOver.sigGO1=hgCondOver.sigGO[1]
    for (hgCondOver.sigGO1 in hgCondOver.sigGO) {
        #hgCondOver.sigGO1= hgCondOver.sigGO[1]
        cat(hgCondOver.sigGO1, "\n")
        hgCondOver1.geneid=geneIdsByCategory(KEGG_hgOver_allGoTerm, catids = hgCondOver.sigGO1)[[1]]
        hgCondOver1.symbol=unlist(mget(hgCondOver1.geneid, org.Mm.egSYMBOL))
        names(hgCondOver1.symbol)=NULL

        if (length(hgCondOver1.symbol)==0){
            #catch case where no geneids found
            hgCondOver.symbol.all=rbind(hgCondOver.symbol.all, "no geneids found")
        }else{
            hgCondOver.symbol.all=rbind(hgCondOver.symbol.all, paste(hgCondOver1.symbol, collapse=", "))
        }       
    }
    
    hgCondOver.sum.all=cbind(KEGG_hgOver.sum_adpv, hgCondOver.symbol.all)
    names(hgCondOver.sum.all)=c(names(KEGG_hgOver.sum_adpv), "Genes")
    summary.file=file.path("/Users/guest1/xin/", "hg_CondOver_sum_all_adpv_KEGG_naiveCocaine_onlycocaine.txt")
    write.table(hgCondOver.sum.all, file=summary.file, row.names=F, col.names=T, sep="\t", quote=F);
    summarygo.file=file.path("/Users/guest1/xin/", "hgCondOver_sum_term_adpv_KEGG_naiveCocaine_onlycocaine.txt")
    write.table(KEGG_hgOver.sum_adpv$Term, file=summarygo.file, row.names=F, col.names=F, sep="\t", quote=F);

    #rm(KEGG_hgOver.sum_adpv, hgCondOver.sum.all)
}