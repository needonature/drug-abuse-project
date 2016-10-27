#source("findDifGenes.R")

data.saline=read.delim(file="/Users/guest1/xin/sig.peakallmet.txt",sep='\t',header=TRUE,stringsAsFactors=FALSE)
data.naive=read.delim(file="/Users/guest1/xin/sig.peakalln.txt",sep='\t',header=TRUE,stringsAsFactors=FALSE)

name.saline.all=name.naive.all=NULL
for (i in 1:nrow(data.saline)){
	if (data.saline[i,1]!=0){
name.saline.all=rbind(name.saline.all,i)
}
}

for (i in 1:nrow(data.naive)){
if(data.naive[i,1]!=0){
name.naive.all=rbind(name.naive.all,i)
}
}
only.naive=setdiff(name.naive.all,name.saline.all)
only.saline=setdiff(name.saline.all,name.naive.all)
inter.groups=intersect(name.naive.all,name.saline.all)

#writting files 
write.table(only.naive,file="/Users/guest1/xin/naiveMethamphetamine_onlynaive.txt",row.name=FALSE,col.names=FALSE)
write.table(only.saline,file="/Users/guest1/xin/naiveMethamphetamine_onlymethamphetamine.txt",row.name=FALSE,col.names=FALSE)
write.table(inter.groups,file="/Users/guest1/xin/naiveMethamphetamine_intergroups.txt",row.name=FALSE,col.names=FALSE)


####################
data.saline=read.delim(file="/Users/guest1/xin/sig.peakallmet.txt",sep='\t',header=TRUE,stringsAsFactors=FALSE)
data.naive=read.delim(file="/Users/guest1/xin/sig.peakalls.txt",sep='\t',header=TRUE,stringsAsFactors=FALSE)

name.saline.all=name.naive.all=NULL
for (i in 1:nrow(data.saline)){
	if (data.saline[i,1]!=0){
name.saline.all=rbind(name.saline.all,i)
}
}

for (i in 1:nrow(data.naive)){
if(data.naive[i,1]!=0){
name.naive.all=rbind(name.naive.all,i)
}
}
only.naive=setdiff(name.naive.all,name.saline.all)
only.saline=setdiff(name.saline.all,name.naive.all)
inter.groups=intersect(name.naive.all,name.saline.all)

#writting files 
write.table(only.naive,file="/Users/guest1/xin/salineMethamphetamine_onlysaline.txt",row.name=FALSE,col.names=FALSE)
write.table(only.saline,file="/Users/guest1/xin/salineMethamphetamine_onlymethamphetamine.txt",row.name=FALSE,col.names=FALSE)
write.table(inter.groups,file="/Users/guest1/xin/salineMethamphetamine_intergroups.txt",row.name=FALSE,col.names=FALSE)
