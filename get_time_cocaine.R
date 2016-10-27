source("library.R")
library(locfit)

#read data
proc.dir="/Users/guest1/xin"
file.data=file.path(proc.dir, "drugabuse_GDS3703.txt")
data.all=read.delim(file.data,stringsAsFactors=F)
data=data.all[, -1]
col.data=colnames(data)

# process original file #
lines.tmp=scan(file="/Users/guest1/xin/GDS3703.soft",sep='\n',what=character(),skip=84,nlines=200)
#rename the column
col.data1=NULL
for (i in 1:ncol(data)){
    array.sample.line=grep(col.data[i], lines.tmp, value=T)
    T=strsplit(lines.tmp[array.sample.line],":")[[1]]
    t=strsplit(T[2],";")[[1]][1]
    a=strsplit(t," ")[[1]][2]
    col.data1=c(col.data1,a)
}
colnames(data)=col.data1

#data.cocaine
data.cocaine=cbind(data[,97:99],data[,1:12])

outfile=cbind(data.all[,1],data.cocaine)
out.file1=gsub(".txt", "cocaine.txt", file.data)
write.table(outfile,out.file1,row.names=FALSE,col.names=TRUE)

#rename data.cocaine
colnames(data.cocaine)<-c("0hr","0hr","0hr","1hr","1hr","1hr","2hr","2hr","2hr","4hr","4hr","4hr","8hr","8hr","8hr")
hours<-rep(NA,ncol(data.cocaine))
hours.test<-c("0hr","0hr","0hr","1hr","1hr","1hr","2hr","2hr","2hr","4hr","4hr","4hr","8hr","8hr","8hr")

for (hour.test in hours.test) {
 hours[grep(hour.test, colnames(data.cocaine))]<-hour.test

}

hours<-as.numeric(sub("hr", "", hours))
nhours=length(hours.test)

#fit my data
idx.peak.all=sig.peak.all=de.peak.all=fold.peak.all=NULL
hours.peak.all=hours.sigpeak.all=NULL

# #scale data: matrix scale:
# data.matrix_scale.cocaine=(data.cocaine-mean(as.matrix(data.cocaine)))/sd(as.matrix(data.cocaine))

# ##use different scale to find time value
# time.matrix_scale.cocaine <- matrix(NA, nrow=nrow(data.cocaine), ncol=5)
# for (i.row in 1:nrow(data.matrix_scale.cocaine)) {
#     # for (i.row in 1:10) {
#     #i.row=1##############
#     if (ceiling(i.row/100)==floor(i.row/100)) cat(i.row, "\n")
#     data.1=data.matrix_scale.cocaine[i.row,]
#     data.my=data.frame("data.1"=as.numeric(data.1), "hours"=hours)    
#     # confidence interval #
#     # plot 95% confidence intervals, with local variance estimate.
#     # function: critfit
#     ## Find bandwidth (nearest neighbor component) using GCV ##
    
#     gcv.my.nn=NULL;
#     part=seq(0.5, 1, by=0.05)
#     for (i in part) {
#         gcv.my.nn=c(gcv.my.nn, gcv(locfit(data.1~lp(hours, deg=2, nn=i),data=data.my))["gcv"])
#     }
# imin.nn=part[which.min(gcv.my.nn)]

# # compute and plot 95% confidence intervals, with m=nhours Bonferonni adjustment
# fit.datamy.optnn=locfit(data.1~lp(hours, deg=2, nn=imin.nn), data=data.my)
# cov.datamy=1-(0.05/nhours) #Bonferonni adjustment #?????????????
# crit(fit.datamy.optnn) <- crit(fit.datamy.optnn,cov=cov.datamy)# change confidence interval

# ## find value at time: 0,1,2,4,8
# time.matrix_scale.cocaine[i.row,]=predict(fit.datamy.optnn, c(0,1,2,4,8)) #da1=preplot(fit.datamy.optnn, where="fitp", band="local",get.data=TRUE)
# }


# #scale data: column scale:
# library(matrixStats)
# data.column_scale.cocaine=(data.cocaine-colMeans(as.matrix(data.cocaine)))/colSds(as.matrix(data.cocaine))

# ##use different scale to find time value
# time.column_scale.cocaine <- matrix(NA, nrow=nrow(data.cocaine), ncol=5)
# for (i.row in 1:nrow(data.column_scale.cocaine)) {
#     # for (i.row in 1:10) {
#     #i.row=1##############
#     if (ceiling(i.row/100)==floor(i.row/100)) cat(i.row, "\n")
#     data.1=data.column_scale.cocaine[i.row,]
#     data.my=data.frame("data.1"=as.numeric(data.1), "hours"=hours)
#     gcv.my.nn=NULL;
#     part=seq(0.5, 1, by=0.05)
#     for (i in part) {
#         gcv.my.nn=c(gcv.my.nn, gcv(locfit(data.1~lp(hours, deg=2, nn=i),data=data.my))["gcv"])
#     }
# imin.nn=part[which.min(gcv.my.nn)]

# # compute and plot 95% confidence intervals, with m=nhours Bonferonni adjustment
# fit.datamy.optnn=locfit(data.1~lp(hours, deg=2, nn=imin.nn), data=data.my)
# cov.datamy=1-(0.05/nhours) #Bonferonni adjustment #?????????????
# crit(fit.datamy.optnn) <- crit(fit.datamy.optnn,cov=cov.datamy)# change confidence interval

# ## find value at time: 0,1,2,4,8
# time.column_scale.cocaine[i.row,]=predict(fit.datamy.optnn, c(0,1,2,4,8)) #da1=preplot(fit.datamy.optnn, where="fitp", band="local",get.data=TRUE)
# }


#find time value
time.cocaine <- matrix(NA, nrow=nrow(data.cocaine), ncol=5)
for (i.row in 1:nrow(data.cocaine)) {
    # for (i.row in 1:10) {
    #i.row=1##############
    if (ceiling(i.row/100)==floor(i.row/100)) cat(i.row, "\n")
    data.1=data.cocaine[i.row,]
    data.my=data.frame("data.1"=as.numeric(data.1), "hours"=hours)
    gcv.my.nn=NULL;
    part=seq(0.5, 1, by=0.05)
    for (i in part) {
        gcv.my.nn=c(gcv.my.nn, gcv(locfit(data.1~lp(hours, deg=2, nn=i),data=data.my))["gcv"])
    }
    imin.nn=part[which.min(gcv.my.nn)]

    # compute and plot 95% confidence intervals, with m=nhours Bonferonni adjustment
    fit.datamy.optnn=locfit(data.1~lp(hours, deg=2, nn=imin.nn), data=data.my)
    cov.datamy=1-(0.05/nhours) #Bonferonni adjustment #?????????????
    crit(fit.datamy.optnn) <- crit(fit.datamy.optnn,cov=cov.datamy)# change confidence interval

    ## find value at time: 0,1,2,4,8
    time.cocaine[i.row,]=predict(fit.datamy.optnn, c(0,1,2,4,8)) #da1=preplot(fit.datamy.optnn, where="fitp", band="local",get.data=TRUE)
}

#scale data: row scale:
library(matrixStats)
time.row_scale.cocaine=(time.cocaine-rowMeans(as.matrix(time.cocaine)))/rowSds(as.matrix(time.cocaine))

#write.table(time.matrix_scale.cocaine,file="/Users/guest1/xin/time_matrix_scale_cocaine.txt",row.name=FALSE,col.names=FALSE)
#write.table(time.column_scale.cocaine,file="/Users/guest1/xin/time_column_scale_cocaine.txt",row.name=FALSE,col.names=FALSE)
write.table(time.row_scale.cocaine,file="/Users/guest1/xin/time_row_scale_cocaine.txt",row.name=FALSE,col.names=FALSE)
