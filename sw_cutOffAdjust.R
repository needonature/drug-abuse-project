#source("cutOffAdjust.R")

library(base)

proc.dir="/Users/StephyTse/Documents/DrugAbuse/addictive_drugs_effect/processeddata/NaiveSalineHeroinMorphine/Difgenes/naiveHeroin/heroin"
file.data=file.path(proc.dir, "kegg_Over_sum_all_cutoff1.txt")
data.all=read.delim(file.data, sep="\t")
pvalue=data.all[,2]
adjp=p.adjust(pvalue, method="hochberg")

adjust.data.all=cbind(data.all[,1:2],adjp,data.all[,3:8])

adj.data.all=NULL
for (i in 1:nrow(adjust.data.all)){
    if (adjust.data.all[i,3]<=0.8){
        adj.data.all=rbind(adj.data.all,adjust.data.all[i,])
    }
}

write.table(adjust.data.all, file="/Users/StephyTse/Documents/DrugAbuse/addictive_drugs_effect/processeddata/NaiveSalineHeroinMorphine/Difgenes/naiveHeroin/heroin/adjust_kegg_Over_sum_all.txt",row.name=FALSE, col.names=TRUE, sep="\t",quote=FALSE)

