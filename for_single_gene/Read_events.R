#!/usr/bin/env Rscript-3.4.1
suppressMessages(library(optparse))

#for command line parsing
args = commandArgs(trailingOnly=TRUE)
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="Event file name from Nanopolish", metavar="character"),
        make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$file) | is.null(opt$out)){
  print_help(opt_parser)
  stop("Input file and output names must be supplied.", call.=FALSE)
}

suppressMessages(library(dplyr))
suppressMessages(library(Rcpp))
suppressMessages(library(pracma))
suppressMessages(library(data.table))

#loading c++ script
Rcpp::sourceCpp("./for_r.cpp")

#loading of event files
dat= fread(paste(opt$file))
print("Done loading")
mod=opt$out

#combine events of same strands and positions
dat.com= dat %>% 
		  mutate(count=round(3012*event_length)) %>% 
		  group_by(contig,read_name,position) %>% 
		  summarise(event_stdv=sd_combine(event_stdv,event_level_mean,count),
					event_level_mean=mean_combine(event_level_mean,count),
					count=sum(count),reference_kmer=unique(reference_kmer))  %>% 
		  ungroup() %>%
		  mutate(event_stdv=ifelse(event_stdv==0,0.01,event_stdv))

#saving the results
assign(paste("dat.f.combined.",mod, sep=""),dat.com)
tmp2=(paste("dat.f.combined.",mod, sep=""))
save(list=(tmp2),file=paste("./dat.f.combined.",mod,".RData", sep=""))
rm(list=ls(pattern="dat.f"))

print("script ran successfully")
