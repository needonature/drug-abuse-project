#read file DifGene file
proc.dir="/Users/guest1/xin"
file.1=file.path(proc.dir,"naiveMethamphetamine_intergroups.txt")
file.2=file.path(proc.dir,"naiveMethamphetamine_onlynaive.txt")
file.3=file.path(proc.dir,"naiveMethamphetamine_onlymethamphetamine.txt")

file.inter=read.delim(file.1,sep="\t",header=FALSE,quote="")
file.heroin=read.delim(file.2,sep="\t",header=FALSE,quote="")
file.naive=read.delim(file.3,sep="\t",header=FALSE,quote="")

#read original data
lines.tmp=scan(file="/Users/guest1/xin/GDS3703.soft",sep='\n',what=character(),skip=196)
Info.tmp=scan(file="/Users/guest1/xin/GPL6105.annot",sep='\n',what=character(),skip=28)


def.inter=def.heroin=def.naive=NULL

for (i in 1:nrow(file.inter)){
	 T=strsplit(lines.tmp[file.inter[i,1]],"\t")[[1]][1]
	 array.sample.line=grep(T,Info.tmp,value=T)
	 t=strsplit(Info.tmp[array.sample.line],"\t")[[1]][4]
	 if (t!=""){
	 def.inter=rbind(def.inter,as.numeric(t))
	}
}
colnames(def.inter)=strsplit(Info.tmp[1],"\t")[[1]][4]
def.inter=unique(def.inter)

for (i in 1:nrow(file.heroin)){
	 T=strsplit(lines.tmp[file.heroin[i,1]],"\t")[[1]][1]
	 array.sample.line=grep(T,Info.tmp,value=T)
	 t=strsplit(Info.tmp[array.sample.line],"\t")[[1]][4]
	 if (t!=""){
	 def.heroin=rbind(def.heroin,as.numeric(t))
	}
}
colnames(def.heroin)=strsplit(Info.tmp[1],"\t")[[1]][4]
def.heroin=unique(def.heroin)

for (i in 1:nrow(file.naive)){
	 T=strsplit(lines.tmp[file.naive[i,1]],"\t")[[1]][1]
	 array.sample.line=grep(T,Info.tmp,value=T)
	 t=strsplit(Info.tmp[array.sample.line],"\t")[[1]][4]
	 if (t!=""){
	 def.naive=rbind(def.naive,as.numeric(t))
	}
}
colnames(def.naive)=strsplit(Info.tmp[1],"\t")[[1]][4]
def.naive=unique(def.naive)

# for (i in 1:nrow(file.naive)){
# 	#i=1
# 	 T=strsplit(lines.tmp[file.naive[i,1]],"\t")[[1]][1]
# 	 array.sample.line=grep(T,Info.tmp,value=T)
# 	 t=strsplit(Info.tmp[array.sample.line],"\t")[[1]][4]
# 	 #if (t!=""){
# 	 if (t!=""){
# 	 #def.naive=rbind(def.naive,as.numeric(t))
# 	 def.naive=rbind(def.naive,as.numeric(t))
# 	}
# }

# def.naive=unique(def.naive)
# colnames(def.naive)=strsplit(Info.tmp[1],"\t")[[1]][4]


out.file1=gsub("naiveMethamphetamine_intergroups.txt","naiveMethamphetamine_inter_geneID.txt",file.1)
out.file2=gsub("naiveMethamphetamine_onlynaive.txt","naiveMethamphetamine_onlynaive_geneID.txt",file.2)
out.file3=gsub("naiveMethamphetamine_onlymethamphetamine.txt","naiveMethamphetamine_onlymethamphetamine_geneID.txt",file.3)
#out.total.id=gsub("nh_onlynaive.txt","total.ID.txt",file.3)

write.table(def.inter,out.file1,col.names=FALSE,row.names=FALSE)
write.table(def.heroin,out.file2,col.names=FALSE,row.names=FALSE)
write.table(def.naive,out.file3,col.names=FALSE,row.names=FALSE)
#write.table(id.all,out.total.id,col.names=TRUE,row.names=FALSE)


########################

#read file DifGene file
proc.dir="/Users/guest1/xin"
file.1=file.path(proc.dir,"salineMethamphetamine_intergroups.txt")
file.2=file.path(proc.dir,"salineMethamphetamine_onlysaline.txt")
file.3=file.path(proc.dir,"salineMethamphetamine_onlymethamphetamine.txt")

file.inter=read.delim(file.1,sep="\t",header=FALSE,quote="")
file.heroin=read.delim(file.2,sep="\t",header=FALSE,quote="")
file.naive=read.delim(file.3,sep="\t",header=FALSE,quote="")

#read original data
lines.tmp=scan(file="/Users/guest1/xin/GDS3703.soft",sep='\n',what=character(),skip=196)
Info.tmp=scan(file="/Users/guest1/xin/GPL6105.annot",sep='\n',what=character(),skip=28)


def.inter=def.heroin=def.naive=NULL

for (i in 1:nrow(file.inter)){
	 T=strsplit(lines.tmp[file.inter[i,1]],"\t")[[1]][1]
	 array.sample.line=grep(T,Info.tmp,value=T)
	 t=strsplit(Info.tmp[array.sample.line],"\t")[[1]][4]
	 if (t!=""){
	 def.inter=rbind(def.inter,as.numeric(t))
	}
}
colnames(def.inter)=strsplit(Info.tmp[1],"\t")[[1]][4]
def.inter=unique(def.inter)

for (i in 1:nrow(file.heroin)){
	 T=strsplit(lines.tmp[file.heroin[i,1]],"\t")[[1]][1]
	 array.sample.line=grep(T,Info.tmp,value=T)
	 t=strsplit(Info.tmp[array.sample.line],"\t")[[1]][4]
	 if (t!=""){
	 def.heroin=rbind(def.heroin,as.numeric(t))
	}
}
colnames(def.heroin)=strsplit(Info.tmp[1],"\t")[[1]][4]
def.heroin=unique(def.heroin)


for (i in 1:nrow(file.naive)){
	 T=strsplit(lines.tmp[file.naive[i,1]],"\t")[[1]][1]
	 array.sample.line=grep(T,Info.tmp,value=T)
	 t=strsplit(Info.tmp[array.sample.line],"\t")[[1]][4]
	 if (t!=""){
	 def.naive=rbind(def.naive,as.numeric(t))
	}
}

def.naive=unique(def.naive)
colnames(def.naive)=strsplit(Info.tmp[1],"\t")[[1]][4]

out.file1=gsub("salineMethamphetamine_intergroups.txt","salineMethamphetamine_inter_geneID.txt",file.1)
out.file2=gsub("salineMethamphetamine_onlysaline.txt","salineMethamphetamine_onlysaline_geneID.txt",file.2)
out.file3=gsub("salineMethamphetamine_onlymethamphetamine.txt","salineMethamphetamine_onlymethamphetamine_geneID.txt",file.3)
#out.total.id=gsub("nh_onlynaive.txt","total.ID.txt",file.3)

write.table(def.inter,out.file1,col.names=FALSE,row.names=FALSE)
write.table(def.heroin,out.file2,col.names=FALSE,row.names=FALSE)
write.table(def.naive,out.file3,col.names=FALSE,row.names=FALSE)
#write.table(id.all,out.total.id,col.names=TRUE,row.names=FALSE)
