# plot 95% confidence intervals, with local variance estimate.
critplot<-function(fit, d=days) {
  x11()
  plot(fit,band="local", get.data=T)
  x11()
  plot(locfit(residuals(fit)~days, deg=2), get.data=T)
}

findpeak<-function(data) {
  ndays=length(data)
  idx.peak=NULL
  i.peak.t=i.val.t=0
  i.day=1; 
  i.next=i.day+1
  while (i.next<=ndays) {
  #
    if (data[i.next]<data[i.day]) {
      if (i.peak.t>i.val.t)
	idx.peak=c(idx.peak, i.peak.t)
      i.val.t=i.next
    } else if (data[i.next]>data[i.day]) {
      if (i.val.t>i.peak.t) idx.peak=c(idx.peak, -i.val.t)
      i.peak.t=i.next 
    }
    i.day=i.day+1
    i.next=i.day+1
  }       
  if (is.null(idx.peak)) {
    if (i.peak.t==i.day) {
      idx.peak=ndays
    } else if (i.val.t==i.day) {
      idx.peak=-ndays
    } else if ((i.peak.t==0)&&(i.val.t==0)) {
      idx.peak=0
    }
  }
  return(idx.peak)
}
