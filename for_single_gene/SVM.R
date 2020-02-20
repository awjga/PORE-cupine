#!/usr/bin/env Rscript-3.4.1

suppressMessages(library(optparse))

args = commandArgs(trailingOnly=TRUE)
option_list = list(
  make_option(c("-u", "--unmod"), type="character", default=NULL,
              help="unmodified RData", metavar="character"),
  make_option(c("-m", "--mod"), type="character", default=NULL,
              help="Modified RData", metavar="character"),
	  make_option(c("-l", "--length"), type="integer", default=NULL,
              help="Modified RData", metavar="integer"),
make_option(c("-o", "--output"), type="character", default=NULL,
              help="output name", metavar="character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$unmod) | is.null(opt$output)|is.null(opt$mod) ){
  print_help(opt_parser)
  stop("Input files and output name must be supplied.", call.=FALSE)
}

suppressMessages(library(dplyr))
suppressMessages(library(e1071))

mod.name=load(paste(opt$mod,sep=""))
dat.mod2=get(mod.name)
rm(list=ls(pattern="dat.f"))
dat.mod2[dat.mod2[,"event_stdv"]==0,"event_stdv"]=0.01

if(NROW(unique(dat.mod2$read_name))<10)
	next()

unmod.name=load(opt$unmod)
dat.unmod_t2=get(unmod.name)
rm(list=ls(pattern="dat.f"))

dat.unmod_t2[dat.unmod_t2[,"event_stdv"]==0,"event_stdv"]=0.01

t=opt$length*0.5
dat.unmod_t=dat.unmod_t2 %>% group_by(read_name)%>%
  filter(n()>t)

dat.mod=dat.mod2 %>% group_by(read_name)%>%
  (n()>t)
  pos=max(dat.mod$position)
  

  pos_m=sort(unique(dat.mod$position))


  mod_mat=matrix("x",nrow=NROW(unique(dat.mod$read_name)),ncol=(pos+1))
  colnames(mod_mat)=0:pos
  rownames(mod_mat)=unique(dat.mod$read_name)
  unmod_mat=matrix(0,ncol=1,nrow=(pos+1))
  nt=matrix("N",ncol=1,nrow=(pos+1))

  for(b in pos_m)
  {
  ### loading of positions
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


	if(NROW(tmpun_t)<5)
		next()

    ### training of unmodified 
	svm.model=svm(tmpun_t[1:2],y=NULL,
			   type='one-classification',
			   nu=9e-04,
			   gamma=0.04,
			   kernel="radial")

    ### prediction of modifications
	svm.predtest=predict(svm.model,(tmp1x[2:3]))
    
    mod_mat[rownames(mod_mat) %in% unlist(tmp1x[,"read_name"]),(b+1)]=0
    mod_mat[rownames(mod_mat) %in% unlist(tmp1x[!svm.predtest,"read_name"]),(b+1)]=1
    unmod_mat[(b+1),]=NROW(tmpun_t)
    nt[(b+1),]=substr(kmer,3,3)
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

write.table(dat.output,paste("Output_",gene_name,opt$output,".csv",sep=""),sep = ",",row.names=FALSE)

rm(list=ls(pattern="dat."))
print("done")
