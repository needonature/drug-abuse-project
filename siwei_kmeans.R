
##source(kmeans.R)
#########################################
kmeans.wei<-function(data.loc, ncl.loc, nstart, iter) {
    
    for (i in 1:iter) {
        
        if (i%%10==0) cat(i,	"\n")
        cl.resloc=kmeans(data.loc, ncl.loc, iter.max=100, nstart=nstart)
        sumss.loc=sum(cl.resloc$withinss)
        cl.loc=cl.resloc$cluster
        
        if (i==1) {
            sumss.minloc=sumss.loc
            cl.minloc=cl.loc
        } else {
            if (sumss.loc<sumss.minloc) {
                sumss.minloc=sumss.loc
                cl.minloc=cl.loc
            }
        }
    }
    res.loc=list(cl.minloc, sumss.minloc)
    names(res.loc)=c("cl", "sumss")
    return(res.loc)
}


##################################################
#library(GLAD)

proc.dir="/Users/StephyTse/Documents/DrugAbuse/addictive_drugs_effect/processeddata/NaiveSalineHeroinMorphine/saline"
file.data=file.path(proc.dir, "saline_fit_data_scalebyrow.txt")
data.all=read.delim(file.data,stringsAsFactors=F,sep="\t")
data=data.all[,-1]

cls=2:6
for (n.cl in cls) { ## n.cl=7 or 8 looks good; n.cl=11 breaks down
    # n.cl=3
    
    cat(n.cl, "\n")
    
    cl.res=i.ord=idxcl.ord=data.ord=data.tv=NULL
    
    cl.res=kmeans.wei(data.loc=data, ncl.loc=n.cl, nstart=20, iter=1000)
    idx.cls=cl.res$cl
    sidx.cls=sort.int(idx.cls, index.return=TRUE)
    sidx.x=sidx.cls$x
    sidx.idx=sidx.cls$ix

cluster=NULL
for (i in 1:length(sidx.x)){
    j=sidx.idx[i]
    clustercol=cbind(sidx.x[i],data.all[j,])
    cluster=rbind(cluster,clustercol)
}


  colnames(cluster)=c("cluster","genesymbol","0","1","2","4","8")
  file.idxcl=file.path(proc.dir,
  gsub("scalebyrow", paste("kIs", n.cl, sep=""), basename(file.data)))
  file.idxcl=gsub(".txt", "tv.txt", file.idxcl)
  write.table(cluster, file.idxcl, row.names=F, col.names=TRUE, sep="\t", quote=F);




}