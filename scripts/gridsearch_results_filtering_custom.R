library(tidyverse)
library(Biobase)
library(QDNAseqmod)
suppressWarnings(library(doMC))
suppressWarnings(library(foreach))

## Added by PS
args = commandArgs(trailingOnly=TRUE)

metafile <- snakemake@params[["meta"]]
metadata <- read.table(file = metafile,header=T,sep="\t", quote="", fill=FALSE)
bin <- as.numeric(snakemake@params[["bin"]])
out_dir <- snakemake@params[["outdir"]]
project <- snakemake@params[["project"]]
purity_cutoff <- as.numeric(snakemake@params[["purity_cutoff"]])
cores <- as.numeric(snakemake@threads)

# results_path <- "/Users/jbreynier/Desktop/Research/Rabadan_Lab/SU2C_organoid/abs_cnv_analysis/absolutecn_output_redo_24_04_17_allsmooth/sWGS_fitting/su2c_organoids_30kb/absolute_PRE_down_sampling/"

# metafile <- "/Users/jbreynier/Desktop/Research/Rabadan_Lab/SU2C_organoid/swgs-absolutecn-mm10/input_data/organoid_sample_sheet_24_04_19_exppurity.tsv"
# metadata <- read.table(file = metafile,header=T,sep="\t", quote="", fill=FALSE)
# bin <- 30
# out_dir <- "/Users/jbreynier/Desktop/Research/Rabadan_Lab/SU2C_organoid/abs_cnv_analysis/absolutecn_output_redo_24_04_17_allsmooth/"
# project <- "su2c_organoids"
# purity_cutoff <- 0.15
# cores <- 2


registerDoMC(cores)

# read in relative CN data
# collapse rds files function
rds.filename <- snakemake@input[["rds"]]

rds.list <- lapply(rds.filename,FUN=function(x){readRDS(x)})

collapse_rds <- function(rds.list){
  comb <- rds.list[[1]][[1]]
  if(length(rds.list) > 1){
    for(i in 2:length(rds.list)){
      add <- rds.list[[i]][[1]]
      comb <- combine(comb,add)
    }
    rds.obj <- comb
  }
  return(rds.obj)
}
# Combine and load rds objects
relative_smoothed <- collapse_rds(rds.list)
saveRDS(relative_smoothed,paste0(out_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_PRE_down_sampling/",project,"_",bin,"kb_relSmoothedCN.rds"))

# relative_smoothed <- readRDS("/Users/jbreynier/Desktop/Research/Rabadan_Lab/SU2C_organoid/abs_cnv_analysis/absolutecn_output_redo_24_04_17_allsmooth/sWGS_fitting/su2c_organoids_30kb/absolute_PRE_down_sampling/su2c_organoids_30kb_relSmoothedCN.rds")

filelist <- snakemake@input[["cl"]]
# filelist <- c("/Users/jbreynier/Desktop/Research/Rabadan_Lab/SU2C_organoid/abs_cnv_analysis/absolutecn_output_redo_24_04_17_allsmooth/sWGS_fitting/su2c_organoids_30kb/absolute_PRE_down_sampling/clonality_results_custom/su2c_organoids_b0_Control_clonality.tsv")
clonality <- do.call(rbind,
			lapply(filelist,FUN = function(x){
				n <- gsub(pattern="_clonality.tsv",rep="",x=x)
        prefix <- paste0(out_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_PRE_down_sampling/clonality_results/",project,"_")
        n <- gsub(pattern=prefix,rep="",x=n)
				tab <- read.table(x,sep="\t",skip=1)
				tab <- cbind(rep(n,times=nrow(tab)),tab)
				return(tab)
			}))
colnames(clonality) <- c("SAMPLE_ID","ploidy","purity","clonality","downsample_depth","powered")

clonality <- left_join(clonality,metadata,by="SAMPLE_ID") %>%
                select(SAMPLE_ID,PATIENT_ID,ploidy,purity,clonality,downsample_depth,powered,exp_p,smooth)
                
## Added by PS
depthtocn<-function(x,purity,seqdepth) #converts readdepth to copy number given purity and single copy depth
{
  (x/seqdepth-2*(1-purity))/purity
}

#Top 10 when ranking by clonality and TP53
#Changed to top 20?

filtered_results <- clonality %>%
  select(SAMPLE_ID, PATIENT_ID, everything()) %>%
  filter(powered ==1) %>%
  group_by(SAMPLE_ID, ploidy) %>%
  mutate(rank_clonality = min_rank(clonality)) %>% #rank clonality within a unique ploidy state 
  filter(rank_clonality ==1) %>% #select ploidy with the lowest clonality within a unique ploidy state 
  group_by(SAMPLE_ID) %>%
  top_n(-20, wt = clonality) %>% # select top 20 ploidy states with the lowest clonality values
  mutate(rank_clonality = min_rank(clonality)) %>% # rank by clonality within a sample across ploidies in top 10
 # retain samples without TP53 mutations and where expected and observed TP53freq <=0.15
  # filter(is.na(TP53freq) | near(expected_TP53_AF,TP53freq, tol = af_cutoff )) %>% 
  filter(is.na(exp_p) | near(purity, exp_p, tol = purity_cutoff )) %>%  
  arrange(PATIENT_ID, SAMPLE_ID)

#Further limit the results by selecting the ploidy states with the lowest clonality values where multiple similar solutions are present.
#Threshold of 0.3 used to select different states

pruned_results <- filtered_results %>%
  arrange(SAMPLE_ID, ploidy) %>%
  group_by(SAMPLE_ID) %>%
  mutate(pl_diff = abs(ploidy-lag(ploidy)), pu_diff = abs(purity-lag(purity))) %>%
  # mutate(new = row_number() == 1 | pl_diff > 0.3) %>%
  mutate(new = row_number() == 1 | pl_diff > 0.3 | pu_diff > 0.3) %>%
  mutate(new_state = cumsum(new)) %>%
  group_by(SAMPLE_ID, new_state) %>%
  filter(rank_clonality == min(rank_clonality))

pruned_results$use <- rep(NA,times=nrow(pruned_results))
pruned_results$notes <- rep(NA,times=nrow(pruned_results))

write.table(filtered_results,
  paste0(out_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_PRE_down_sampling/",project,"_filtered_results.tsv"),
  sep="\t",col.names=T,row.names=F,quote=F)

write.table(pruned_results,
  paste0(out_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_PRE_down_sampling/",project,"_fit_QC_predownsample.tsv"),
  sep="\t",col.names=T,row.names=F,quote=F)

## ADDED by PS - adding output folder for results
if(!dir.exists(paste0(out_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_PRE_down_sampling/plots"))){
	dir.create(paste0(out_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_PRE_down_sampling/plots"))
}

#relative_smoothed
#Plot absolute CN fits for assessment
foreach(i=unique(pruned_results$SAMPLE_ID)) %dopar% {
  dat <-  pruned_results %>%
    filter(SAMPLE_ID == i) %>%
    arrange(ploidy)
    #arrange(rank_clonality)
  x <- relative_smoothed[, i]
  cn <- assayDataElement(x,"copynumber")
  seg <- assayDataElement(x,"segmented")
  rel_ploidy <- mean(cn,na.rm=T)
  ll <- nrow(dat)
  png(paste0(out_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_PRE_down_sampling/plots/", i, ".png"), w= 450*ll, h = 350)
  par(mfrow = c(1,ll))
  for(n in 1:nrow(dat)){
    
    ploidy <- dat[n,]$ploidy
    purity <- dat[n,]$purity
    cellploidy <- ploidy*purity+2*(1-purity)
    seqdepth <- rel_ploidy/cellploidy
    
    # expTP53 <- round(dat[n,]$expected_TP53_AF, 2)
    # TP53 <- dat[n,]$TP53freq
    exp_purity <- dat[n,]$exp_p
    
    #convert to abs
    
    pData(x)$ploidy <- ploidy
    pData(x)$purity <- purity
    
    temp <- x
    abs_cn <- depthtocn(cn,purity,seqdepth)
    abs_seg <- depthtocn(seg,purity,seqdepth)
    assayDataElement(temp,"copynumber") <- abs_cn
    assayDataElement(temp,"segmented") <- abs_seg
    
    #tmp_abs <- convert_rd_to_cn(x)
 # plot   
    if(ploidy>5){
      yrange=15
    }else
    {
      yrange=10
    }
  plot(temp,doCalls=FALSE,showSD=TRUE,logTransform=FALSE,ylim=c(0,yrange),ylab="Absolute tumour CN",
         main=paste(i,
                    " exp_p=", round(exp_purity,2),
                    " p=",round(purity,2),
                    " pl=",round(ploidy,2),
                    sep=""),cex.main=0.8)
  abline(h=1:9,col = "blue")
  
  }
  dev.off()
}
