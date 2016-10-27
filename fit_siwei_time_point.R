library(locfit)
source("library.R")

#read data
proc.dir="/Users/guest1/siwei/heroin"
file.data=file.path(proc.dir, "heroine.txt")
data.all=read.delim(file.data,stringsAsFactors=F,sep=" ")
lines.tmp=scan(file="/Users/guest1/xin/GDS3703.soft",sep='\n',what=character(),skip=196)



data.heroin=data.all[,-1]

colnames(data.heroin)<-c("0hr","0hr","0hr","1hr","1hr","1hr","2hr","2hr","2hr","4hr","4hr","4hr","8hr","8hr","8hr")
hours<-rep(NA,ncol(data.heroin))
hours.test<-c("0hr","0hr","0hr","1hr","1hr","1hr","2hr","2hr","2hr","4hr","4hr","4hr","8hr","8hr","8hr")

for (hour.test in hours.test) {
 hours[grep(hour.test, colnames(data.heroin))]<-hour.test

}

hours<-as.numeric(sub("hr", "", hours))
nhours=length(hours.test)


idx.peak.all=sig.peak.all=de.peak.all=fold.peak.all=NULL
hours.peak.all=hours.sigpeak.all=NULL


fit.data.all=fit.data.scalerow.all=fit.data.only=NULL
for (i.row in 1:nrow(data.heroin)) {
    # for (i.row in 1:10) {
    # i.row=3
    if (ceiling(i.row/100)==floor(i.row/100)) cat(i.row, "\n")
    data.1=data.heroin[i.row,]
    data.my=data.frame("data.1"=as.numeric(data.1), "hours"=hours)
    
    
    # confidence interval #
    # plot 95% confidence intervals, with local variance estimate.
    # function: critfit
    
    
    ## Find bandwidth (nearest neighbor component) using GCV ##
    
    gcv.my.nn=NULL;
    part=seq(0.5, 1, by=0.05)
    for (i in part) {
        gcv.my.nn=c(gcv.my.nn, gcv(locfit(data.1~lp(hours, deg=2, nn=i),data=data.my))["gcv"])
    }

# x11()
# plot(part, gcv.my.nn)
imin.nn=part[which.min(gcv.my.nn)]


fit.datamy.optnn=locfit(data.1~lp(hours, deg=2, nn=imin.nn), data=data.my)
cov.datamy=1-(0.05/nhours) #Bonferonni adjustment #?????????????
crit(fit.datamy.optnn) <- crit(fit.datamy.optnn,cov=cov.datamy)

m.data=preplot(fit.datamy.optnn, band="local",get.data=TRUE) ##cannot find preplot.locafit????

fit.data=m.data$fit
fit.se=m.data$se.fit
fit.hours=m.data$xev[[1]]
crit.const=m.data$critval$crit.val
ci.top=fit.data+crit.const * fit.se
ci.low=fit.data-crit.const * fit.se

##
## The negative and positive value of a peak in idx.peaks
## is determined by the value of this peak w.r.t. that of
## its predecessor
## If
##   the value of the peak > predecessor: positive
##   the value of the peak < predecessor: negative
##

idx.peak=sig.peak=de.peak=NULL

idx.peak=findpeak(fit.data)
idx.peak.all=c(idx.peak.all, list(idx.peak))


## Identify significant peaks ##

n.peak=length(idx.peak)
sig.peak=NULL
if ((n.peak==1)&&(idx.peak==0)) {
    ### EXIT!!! No DE###
    sig.peak=0
} else {
    
    while (is.null(sig.peak)&&(n.peak>0)) {
        i.peak=1
        idx.peak1=idx.peak[i.peak]
        ci.low.idxpeak1=ci.low[abs(idx.peak1)]
        ci.top.idxpeak1=ci.top[abs(idx.peak1)]
        fit.peak=fit.data[abs(idx.peak1)]
        if ((idx.peak1<0)&&(ci.top.idxpeak1<ci.low[1])&&((fit.data[1]-fit.peak)>log2(1.2))) {
            sig.peak=c(sig.peak, idx.peak1)
            # break
        } else if ((idx.peak1>0)&&(ci.low.idxpeak1>ci.top[1])&&((fit.peak-fit.data[1])>log2(1.2))) {
            sig.peak=c(sig.peak, idx.peak1)
            # break
        } else if (n.peak>1) {
            idx.peak=idx.peak[-1]
            n.peak=length(idx.peak)
        } else if (n.peak==1) {
            sig.peak=0
        }
    } # end while
    
    if (n.peak>1) {
        for (i.peak in 2:n.peak) {
            #
            idx.peak1=idx.peak[i.peak]
            ci.low.idxpeak1=ci.low[abs(idx.peak1)]
            ci.top.idxpeak1=ci.top[abs(idx.peak1)]
            fit.peak=fit.data[abs(idx.peak1)]
            idx.peak.pre=idx.peak[(i.peak-1)]
            if ((idx.peak1<0)&&(ci.low.idxpeak1<ci.top[abs(idx.peak.pre)])&&((fit.data[idx.peak.pre]-fit.peak)>=log2(1.2))) {
                sig.peak=c(sig.peak, idx.peak1)
            }else if ((idx.peak1>0)&&(ci.low.idxpeak1>ci.top[abs(idx.peak.pre)])&&((fit.peak-fit.data[idx.peak.pre])>=log2(1.2))) {
                sig.peak=c(sig.peak, idx.peak1)
            }
    
            
        } # end for i.peak in 2:n.peak
    } # end if (n.peak>1)
}

if (sig.peak!=0){
fit.predictdata=predict(fit.datamy.optnn,c(0,1,2,4,8))
fit.data.scalebyrow=scale(fit.predictdata)

T=strsplit(lines.tmp[i.row],"\t")[[1]]
t=paste("genes",i.row,T[2])
fit.data.line=cbind(T[1],t,t(fit.predictdata))
fit.data.scaledrow.line=cbind(T[1],t,t(fit.data.scalebyrow))
fit.data.scalerow.all=rbind(fit.data.scalerow.all,fit.data.scaledrow.line)
fit.data.all=rbind(fit.data.all, fit.data.line)
fit.data.only=rbind(fit.data.only,t(fit.predictdata))

}
}
fit.scalecol.all=data.std=NULL
for (i in 1:ncol(fit.data.only)){
	fit.data.col=scale(fit.data.only[,i])
	fit.scalecol.all=cbind(fit.scalecol.all,fit.data.col)
}

fit.scalecolumn.all=cbind(fit.data.all[,1:2],fit.scalecol.all)

average=mean(as.matrix(fit.data.only))
std=sd(as.matrix(fit.data.only))

for (i in 1:nrow(fit.data.only)){
    if (i%%1000==0)cat(i, "\n")
    data.column.all=NULL
    for (j in 1:ncol(fit.data.only)){
        data.column=(fit.data.only[i,j]-average)/std
        data.column.all=c(data.column.all,data.column)
    }
    data.std=rbind(data.std,data.column.all)
}

data.std.all=cbind(fit.data.all[,1:2],data.std)

colnames(fit.data.all)=c("ID","genesymbol","0h","1h","2","4","8")
colnames(fit.data.scalerow.all)=c("ID","genesymbol","0h","1h","2","4","8")
colnames(fit.scalecolumn.all)=c("ID","genesymbol","0h","1h","2","4","8")
colnames(data.std.all)=c("ID","genesymbol","0h","1h","2","4","8")

write.table(fit.data.all, file="/Users/guest1/siwei/heroin/heroin_fit_data.txt",sep="\t",row.names=FALSE,col.names=TRUE)
write.table(fit.data.scalerow.all, file="/Users/guest1/siwei/heroin/heroin_fit_data_scalebyrow.txt",sep="\t",row.names=FALSE,col.names=TRUE)
write.table(fit.scalecolumn.all, file="/Users/guest1/siwei/heroin/heroin_fit_data_scalebycol.txt",sep="\t",row.names=FALSE,col.names=TRUE)
write.table(data.std.all, file="/Users/guest1/siwei/heroin/heroin_fit_data_scalebymatrix.txt",sep="\t",row.names=FALSE,col.names=TRUE)


