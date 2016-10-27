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


# Get the entrez gene identifiers that are mapped to a GO ID
x <- org.Mm.egGO
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes]) 

#find our geneID in geneID_all
inds <- which(names(xx) %in% geneID_all) 
gene_all<-names(xx[inds])

#use annot file to find geneID of clustered gene
proc.dir="/Users/guest1/xin"
file.data=file.path(proc.dir, "GPL6105_revised.annot")
annotation=read.delim(file.data,stringsAsFactors=F)

for (n in 2:15) {

    proc.dir="/Users/guest1/xin"
    file.data=file.path(proc.dir, paste0("cocaine_",paste("kIs", n, sep=""),"tv.txt"))
    clusterN=read.delim(file.data,stringsAsFactors=F)

    for (c in 1:n) {
        clusterC<-clusterN[which(clusterN$cluster == c),]
        
        #find our clusterN_gene
        l<-strsplit(clusterC$genesymbol, "_")
        clusterN_gene<-NA
        for(ii in 1:length(l)){
            clusterN_gene<-c(clusterN_gene,l[[ii]][1])
        }
        clusterN_gene<-clusterN_gene[-1]
        clusterN_gene<-unique(clusterN_gene)########

        #find our clusterN_geneID from annotation
        clusterN_geneID<-annotation[as.numeric(clusterN_gene),]$Gene.ID
        clusterN_geneID<-clusterN_geneID[which(clusterN_geneID!="")]
        inds <- which(gene_all %in% clusterN_geneID) 
        clusterN_geneID<-gene_all[inds]

        #conditional GOHyperGeo
        hgCutoff <- 0.01  #when we set 1, we get all the GO term's p-value
        params <- new("GOHyperGParams",
                      geneIds=clusterC,
                      universeGeneIds=gene_all,
                      annotation="org.Mm.eg.db",
                      ontology="BP",
                      pvalueCutoff=hgCutoff,
                      conditional=TRUE,
                      testDirection="over")
        GO_hgCondOver <- hyperGTest(params)
        GO_hgCondOver
        GO_hgCondOver.sum=summary(GO_hgCondOver)

        #for adjusting p values, first find all go term
        pvalueCutoff(paramsCond) <- 1
        GO_hgCondOver_allGoTerm <- hyperGTest(paramsCond)
        GO_hgCondOver.sum_adpv<-summary(GO_hgCondOver_allGoTerm)
        summary.file=file.path("/Users/guest1/xin/", paste0("hgCondOver_sum_all_allGOterm_GO_cocaine_",paste(n, c, sep="_"),".txt"))
        write.table(GO_hgCondOver.sum_adpv, file=summary.file, row.names=F, col.names=T, sep="\t", quote=F);
        summarygo.file=file.path("/Users/guest1/xin/", paste0("hgCondOver_sum_term_allGOterm_GO_cocaine_",paste(n, c, sep="_"),".txt"))
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
                #cat(hgCondOver.sigGO1, "\n")
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
            summary.file=file.path("/Users/guest1/xin/", paste0("hgCondOver_sum_all_adpv_GO_cocaine_",paste(n, c, sep="_"),".txt"))
            write.table(hgCondOver.sum.all, file=summary.file, row.names=F, col.names=T, sep="\t", quote=F);
            summarygo.file=file.path("/Users/guest1/xin/", paste0("hgCondOver_sum_term_adpv_GO_cocaine_",paste(n, c, sep="_"),".txt"))
            write.table(GO_hgCondOver.sum_adpv$Term, file=summarygo.file, row.names=F, col.names=F, sep="\t", quote=F);

            #rm(GO_hgCondOver.sum_adpv, hgCondOver.sum.all)
        }
    }
}