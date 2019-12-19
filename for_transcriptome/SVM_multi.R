#!/usr/bin/env Rscript-3.4.1

suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(e1071))
suppressMessages(library(data.table))


args = commandArgs(trailingOnly=TRUE)
option_list = list(
  make_option(c("-u", "--unmod"), type="character", default=NULL,
              help="unmodified RData", metavar="character"),
  make_option(c("-m", "--mod"), type="character", default=NULL,
              help="Modified RData", metavar="character"),
	make_option(c("-s", "--split"), type="integer", default=NULL,
		help="Event file name from Nanopolish", metavar="integer"),
	make_option(c("-p", "--portion"), type="integer", default=NULL,
		help="Number of cores allocated", metavar="integer"),
	make_option(c("-f", "--folder"), type="character", default=NULL,
		help="Folder name", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$split) | is.null(opt$portion)| is.null(opt$folder))
{
	print_help(opt_parser)
	stop("Folder name must be supplied", call.=FALSE)
}

path_listmod=list.files(opt$mod,pattern=".RData$")
path_listunmod=list.files(opt$unmod,pattern=".RData$")


tmpmod=strsplit(path_listmod,split="combined")
tmpun=strsplit(path_listunmod,split="combined")

tmpmod1=do.call(rbind.data.frame, tmpmod)[,2]
tmpun1=do.call(rbind.data.frame, tmpun)[,2]

tmp_inter=tmpun1[tmpun1 %in% tmpmod1]

tmp_split=round(length(tmp_inter)/opt$split,0)
tmpstart=1+(tmp_split*(opt$portion-1))
tmp_end=tmp_split*(opt$portion)

if(opt$split==opt$portion)
	tmp_end=length(tmp_inter)

dir.create(paste(opt$folder,sep="_"))

for (i in tmpstart:tmp_end){

	gene=as.character(tmp_inter[i])
	mod.name=load(paste(opt$mod,"/",path_listmod[tmpmod1==gene],sep=""))
	dat.mod=get(mod.name)
	rm(list=ls(pattern="dat.f"))

	unmod.name=load(paste(opt$unmod,"/",path_listunmod[tmpun1==gene],sep=""))
	dat.unmod_t=get(unmod.name)
	rm(list=ls(pattern="dat.f"))

	t=dat.unmod_t %>% group_by(read_name) %>% summarise(n=n())
	dat.unmod_t=dat.unmod_t %>% group_by(read_name)%>%
		filter(n()>=(quantile(t$n,0.5)))

	dat.mod=dat.mod %>% group_by(read_name)%>%
		filter(n()>=(quantile(t$n,0.5)))

	gene_name=unique(dat.mod$contig)

	### filter by length
	if(NROW(unique(dat.mod$read_name))<10)
		next()

	pos=max(dat.mod$position)
	pos_m=sort(unique(dat.mod$position))
	mod_mat=matrix("x",nrow=NROW(unique(dat.mod$read_name)),ncol=(pos+1))
	colnames(mod_mat)=0:pos
	rownames(mod_mat)=unique(dat.mod$read_name)
	nt=matrix("N",ncol=1,nrow=(pos+1))
	unmod_mat=matrix(0,ncol=1,nrow=(pos+1))

	for(b in pos_m)
	{
		tmp1x= dat.mod%>%
			filter(position == b) %>%
			select(read_name,event_level_mean,event_stdv,count,reference_kmer)%>%
			mutate(event_level_mean=(event_level_mean),event_stdv=log(event_stdv),count=log(count))

		kmer=first(tmp1x$reference_kmer)
		
		tmpun_t= dat.unmod_t %>%
			ungroup() %>%
			filter(position == b) %>%
			select(event_level_mean,event_stdv,count) %>%
			mutate(event_level_mean=(event_level_mean),event_stdv=log(event_stdv),count=log(count))

		if(NROW(tmpun_t)<20)
			next()

		svm.model=svm(tmpun_t[1:2],y=NULL,
			type='one-classification',
			nu=9e-4,
			gamma=0.04,
			kernel="radial")

		 svm.predtest=predict(svm.model,(tmp1x[2:3]))

		mod_mat[rownames(mod_mat) %in% unlist(tmp1x[,"read_name"]),(b+1)]=0
		mod_mat[rownames(mod_mat) %in% unlist(tmp1x[!svm.predtest,"read_name"]),(b+1)]=1

		nt[(b+1),]=substr(kmer,3,3)
		unmod_mat[(b+1),]=NROW(tmpun_t)

	}
	tmp_mat=mod_mat
	mode(tmp_mat)="numeric"
	mod_strands=as_tibble(apply(tmp_mat, 2, function (x) sum(!is.na(x))))


	gene_name=unique(dat.mod$contig)

	dat.results=as_tibble(colSums(tmp_mat,na.rm=T))
	dat.percentage=as_tibble(dat.results/mod_strands)
	gene_mat=matrix(gene_name,ncol=1,nrow=(pos+1))

	dat.output=bind_cols(as_tibble(gene_mat),as_tibble(3:(pos+3)),as_tibble(nt),dat.percentage,dat.results,mod_strands,as_tibble(unmod_mat))

	colnames(dat.output)=c("Gene","Position","NT","Mod_percentage","Mod_number","Mod_strands","Training_strands")
	if(sum(is.na(dat.percentage))==NROW(dat.percentage))
		next()

	write.table(dat.output,paste(opt$folder,"/Output_",gene_name,".csv",sep=""),sep = ",",row.names=FALSE)
	rm(list=ls(pattern="dat."))
}
print("done")

