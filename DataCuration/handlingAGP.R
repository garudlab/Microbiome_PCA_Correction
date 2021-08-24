mainDir <- "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/"
require(dplyr)
dir1 <- paste0(mainDir,"AGP_max_k5/")
dir2 <- paste0(mainDir,"AGP_complete_otu/")

new_dir1 <- paste0(mainDir,"AGPr_max_k5/")
new_dir2 <- paste0(mainDir,"AGPr_complete_otu/")

kmer_table = read.csv(paste0(dir1,"trim.txt"), sep = "\t", stringsAsFactors = FALSE,header=TRUE,row.names=1)
otu_table = read.csv(paste0(dir2,"trim.txt"),  sep = "\t", stringsAsFactors = FALSE,header=TRUE,row.names=1)
kmer_metadata = read.csv(paste0(dir1,"metadata.txt"), sep = "\t",stringsAsFactors = FALSE,header=TRUE,row.names=1)

otu_metadata = read.csv(paste0(dir2,"metadata.txt"), sep = "\t",stringsAsFactors = FALSE,header=TRUE,row.names=1)

## NOTES OF NEEDS
# 1 filer out repeats
# 2 only feces 
# 3 make sure kmers and outs have same samples

# FINAL FILTER
dim(otu_metadata)
otu_metadata = otu_metadata %>% filter(body_habitat.x == "UBERON:feces",!grepl("BLANK",sample_name),!grepl("Nextera",sample_plate ),
                                        !grepl("resubmission|Nextera", center_project_name))
otu_metadata = otu_metadata %>% arrange(Instrument,desc(collection_month_year))
otu_metadata = otu_metadata[!duplicated(otu_metadata$host_subject_id.x), ]


kmer_metadata = kmer_metadata %>% filter(sample_name %in% row.names(otu_metadata))
kmer_metadata = kmer_metadata %>% filter(body_habitat.x == "UBERON:feces",!grepl("BLANK",sample_name),!grepl("Nextera",sample_plate ),
                                       !grepl("resubmission|Nextera", center_project_name))
kmer_metadata = kmer_metadata %>% arrange(Instrument,desc(collection_month_year))
kmer_metadata = kmer_metadata[!duplicated(kmer_metadata$host_subject_id.x), ]

#sort(table(kmer_metadata$host_subject_id.x),decreasing=TRUE)[1:10]
#sort(table(kmer_metadata$sample_name),decreasing=TRUE)[1:10]
#table(kmer_metadata$bin_antibiotic_last_year,kmer_metadata$Instrument)
dir.create(new_dir1)
dir.create(new_dir2)
write.table(kmer_metadata,paste0(new_dir1,"metadata.txt"),quote = FALSE,sep = "\t")
write.table(otu_metadata,paste0(new_dir2,"metadata.txt"),quote = FALSE,sep = "\t")

#####


## $$$$$ ##### NOTES OF OBSERVATIONS
# 1 19612 samples feces out of 24903 in kmer
# 2 15073 samples out of 19386 in otu
# 3 min(relative abundance)/2 pseudocount
# 4 kmer: AGP_Bloom_Nextera  378 samples  (Sejin Song The results, spanning 15 individuals and over 1,200 samples)
# 5 kmer: BLANK 2188 samples grepl("BLANK",otu_metadata$sample_name)
# 6 otu: BLANK 1761 samples grepl("BLANK",otu_metadata$sample_name)
# 7 otu: 1388 individuals with multiple sequencings,  5513 samples
# 8 otu: only 47/928 samples in Instrument HiSeq 2500 in MiSeq
# 9 otu: only 26/ 413 HiSeq 2000 in HiSeq
# 10: finally 12952  samples # actually before removing duplicated check for overlap with kmer
# 
### $$$$ ####

