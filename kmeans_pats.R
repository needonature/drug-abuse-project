# source("kmeans_pats.R")

############################
## do Kmeans for patients ##
############################

### function ###

kmeans.wei<-function(data.loc, ncl.loc, nstart, iter=100) {
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


# ################

# rm(list=ls())

# library(GLAD)

# source("/Users/weiwu2/Documents/Sally/SARP3/code/dir_info.R")
# source(libfile1)
# source(libfile2)

# # read SARP3 scale data #

# file.data=file.path(proc.dir, "scalesarp150918sptr1_kmeans5_v0tv.txt")
# data.all=read.delim(file.data, stringsAsFactors = F)
# data=apply(data.all[, 3:ncol(data.all)], 2, as.numeric)

# RBcol <- myPalette(low="blue3", high="red3", k=100)

# # get parameters #

# results=c("results1", "results2", "results3")

# for (result in results) {
# # result="results1"

#  cat("Processing", result, "\n")

#  cl.dir=file.path(kmeans.dir, result)
#  if (!file.exists(cl.dir))
#   dir.create(cl.dir)

#  idx.cls=NULL
#  cls=4:8
#  for (n.cl in cls) { ## n.cl=7 or 8 looks good; n.cl=11 breaks down
#  # n.cl=6

#   cat(n.cl, "\n")

#   cl.res=i.ord=idxcl.ord=data.ord=data.tv=NULL

#   cl.res=kmeans.wei(data.loc=data, ncl.loc=n.cl, nstart=20, iter=1000)

#   idx.cls=cbind(idx.cls, cl.res$cl)

#   file.idxcl=file.path(kmeans.dir, result,
#             gsub("scalesarp150918sptr1", paste("idxsamples_cl", n.cl, sep=""), basename(file.data)))
#   file.idxcl=gsub("tv.txt", ".txt", file.idxcl)
#   write.table(cl.res$cl, file.idxcl, row.names=F, col.names=F, sep="\t", quote=F);

#   i.ord=order(cl.res$cl)
#   idxcl.ord=cl.res$cl[i.ord]
#   data.ord=data[i.ord, ]

#   data.tv=data.all[i.ord, ]
#   file.cl.tv=file.path(kmeans.dir, result,
#             gsub("scalesarp150918sptr1", paste("scaledata_samplescl", n.cl, sep=""), basename(file.data)))
#   write.table(data.tv, file.cl.tv, row.names=F, col.names=T, sep="\t", quote=F);

#   ## plot affinity matrix for clustered data ##

#   affinity.data=as.matrix(dist(data.ord, diag=T, upper=T));

#   file.affi=file.path(kmeans.dir, result,
#             gsub("scalesarp150918sptr1", paste("affinitysamples_cl", n.cl, sep=""), basename(file.data)))
#   file.affi=gsub("tv.txt", ".pdf", file.affi)
#   pdf(file.affi)
#   # x11()
#   image(affinity.data, col=RBcol)
#   title('Image of kmeans-clustered affinity matrix of samples')
#   dev.off()

#  } # end for (n.cl in cls ...

#  idx.cls.all=cbind(data.all[, 1:2], idx.cls)
#  colnames(idx.cls.all)=c(colnames(data.all)[1:2], paste("Kmeans_", cls, "cl", sep=""))
#  file.idxcls.all=file.path(kmeans.dir, result, gsub("scalesarp150918sptr1", "idxsamples_allcls",basename(file.data)))
#  file.idxcls.all=gsub("tv.txt", ".txt", file.idxcls.all)
#  write.table(idx.cls.all, file.idxcls.all, row.names=F, col.names=T, sep="\t", quote=F);

# }