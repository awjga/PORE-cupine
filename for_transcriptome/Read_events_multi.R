#!/usr/bin/env Rscript-3.4.1

suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(Rcpp))
suppressMessages(library(data.table))



args = commandArgs(trailingOnly=TRUE)
option_list = list(
  make_option(c("-s", "--split"), type="integer", default=NULL,
              help="number of tread to run", metavar="integer"),
  make_option(c("-p", "--portion"), type="integer", default=NULL,
              help="portion", metavar="integer")
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="folder input", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output folder name", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$split) | is.null(opt$portion)| is.null(opt$input)| is.null(opt$output)){
  print_help(opt_parser)
  stop("Number of parts must be provided", call.=FALSE)
}


### loading c functions
Rcpp::sourceCpp("./for_r.cpp")

### loading files
path_list=list.files(opt$input, recursive=T,pattern=".tmp$")

### spliting files to number of allocated parts
tmp_split=round(length(path_list)/opt$split,0)

### assigning parts to run
tmpstart=1+(tmp_split*(opt$portion-1))
tmp_end=tmp_split*(opt$portion)
### reading header file 
header=fread(paste(opt$input,"/header",sep=""))
if(opt$split==opt$portion)
tmp_end=length(path_list)

dir.create(opt$output)


for (i in tmpstart:tmp_end)
{
  dat=as_tibble(fread(paste(opt$input,"/",path_list[i],sep="")))
  colnames(dat)=colnames(header)

#combine events of same strands and positions
  tmp_dat=dat %>%  
	  mutate(count=round(3012*event_length)) %>%
	  group_by(read_name,position) %>%
	  summarise(contig=first(contig), event_stdv=sd_combine(event_stdv,event_level_mean,count),
				event_level_mean=mean_combine(event_level_mean,count),
				count=sum(count),reference_kmer=unique(reference_kmer))  %>%
	  ungroup() %>%
	  mutate(event_stdv=ifelse(event_stdv==0,0.01,event_stdv))
  
  gene_name=unique(tmp_dat$contig)
  
  #saving the results
  assign(paste("dat.f.combined.",gene_name, sep=""),tmp_dat)
  tmp2=(paste("dat.f.combined.",gene_name, sep=""))
save(list=(tmp2),file=paste(opt$output,"/dat.f.combined.",gene_name,".RData", sep=""))
rm(list=ls(pattern="dat.f"))
}

print("done")
