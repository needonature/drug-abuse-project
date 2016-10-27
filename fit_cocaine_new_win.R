source("library.R")
library(locfit)

#read data
proc.dir="C:\\Users\\thinkpad\\Desktop\\drug abuse project"
file.data=file.path(proc.dir, "drugabuse_GDS3703.txt")
data.all=read.delim(file.data,stringsAsFactors=F)
data=data.all[, -1]
col.data=colnames(data)

# process original file #
lines.tmp=scan(file="C:\\Users\\thinkpad\\Desktop\\drug abuse project\\GDS3703.soft",sep='\n',what=character(),skip=84,nlines=200)
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

######
#plot pdf
pdf(file="C:\\Users\\thinkpad\\Desktop\\drug abuse project\\peaks.pdf",width = 6.2, height = 8,pointsize = 12, onefile = TRUE)
par(mfrow=c(4,3))

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
for (i.row in 1:nrow(data.cocaine)) {
    # for (i.row in 1:10) {
    # i.row=8
    if (ceiling(i.row/100)==floor(i.row/100)) cat(i.row, "\n")
    data.1=data.cocaine[i.row,]
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
# i.row=2; imin.nn=0.75

# compute and plot 95% confidence intervals, with m=nhours Bonferonni adjustment

fit.datamy.optnn=locfit(data.1~lp(hours, deg=2, nn=imin.nn), data=data.my)
cov.datamy=1-(0.05/nhours) #Bonferonni adjustment #?????????????
crit(fit.datamy.optnn) <- crit(fit.datamy.optnn,cov=cov.datamy)# change confidence interval


# plot confidence intervals
# critplot(fit.datamy.optnn)


## Find peaks ##
## Find bumps and valleys ##

m.data=preplot(fit.datamy.optnn, band="local",get.data=TRUE)
fit.data=m.data$fit
fit.se=m.data$se.fit
fit.hours=m.data$xev[[1]]
crit.const=m.data$critval$crit.val
ci.top=fit.data+crit.const * fit.se
ci.low=fit.data-crit.const * fit.se

##plot the peaks


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
            
            idx.peak.pre=idx.peak[(i.peak-1)]
            if ((idx.peak1<0)&&(ci.low.idxpeak1<ci.top[abs(idx.peak.pre)])&&((fit.data[idx.peak.pre]-fit.peak)>=log2(1.2))) {
                sig.peak=c(sig.peak, idx.peak1)
            }else if ((idx.peak1>0)&&(ci.low.idxpeak1>ci.top[abs(idx.peak.pre)])&&((fit.peak-fit.data[idx.peak.pre])>=log2(1.2))) {
                sig.peak=c(sig.peak, idx.peak1)
            }
            
        } # end for i.peak in 2:n.peak
    } # end if (n.peak>1)
} # else (idx.peak!=0)

sig.peak.all=c(sig.peak.all, list(sig.peak))

# hours when sig.peak appear
hours.sigpeak=fit.hours[abs(sig.peak)]*(sig.peak/abs(sig.peak))
hours.sigpeak.all=c(hours.sigpeak.all, list(hours.sigpeak))

## Identify peaks differentially expressed w.r.t. controls ##
##
## The negative and positive value of a peak in de.peak
## is determined by the value of the peak w.r.t. that of
## the control
## If
##   the value of the peak > control: positive
##   the value of the peak < control: negative
##

nsig.peak=length(sig.peak)
# ci.low[1]: ci.low for control
# ci.top[1]: ci.top for control
de.peak=rep(NA, nsig.peak)

if ((length(sig.peak)==1)&&(sig.peak==0)) {
    de.peak=0
} else {
    for (isig.peak in 1:nsig.peak) {
        # isig.peak=1
        sig.peak1=sig.peak[isig.peak]
        data.sigpeak1=fit.data[abs(sig.peak1)]
        ci.low.sigpeak1=ci.low[abs(sig.peak1)]
        ci.top.sigpeak1=ci.top[abs(sig.peak1)]
        
        if ((ci.low.sigpeak1>ci.top[1])&&((data.sigpeak1-fit.data[1])>log2(1.2))) {
            de.peak[isig.peak]=1
        } else if ((ci.top.sigpeak1<ci.low[1])&&((fit.data[1]-data.sigpeak1)>log2(1.2))) {
            de.peak[isig.peak]=-1
        } else {
            de.peak[isig.peak]=0
        }
    }
}
de.peak.all=c(de.peak.all, list(de.peak))
    
   if (sig.peak[1]!=0){
    plot(m.data,ylab=data.all[i.row,2])
}

} # end for irow.data

dev.off()



source("listIO.R")
writelist(idx.peak.all, filename="idx.peakallc.txt")
writelist(sig.peak.all, filename="sig.peakallc.txt")
writelist(de.peak.all, filename="de.peakallc.txt")


# write fit.hours #
#fit.hours.f=as.numeric(formatC(fit.hours))
#write.table(fit.hours.f, file="/Users/StephyTse/Documents/DrugAbuse/addictive_drugs_effect/processeddata/NSH/fithours.txt", sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)