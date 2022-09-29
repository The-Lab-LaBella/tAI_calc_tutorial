# tAI_calc_tutorial
Description of my method using tAI.R (https://github.com/mariodosreis/tai)

I have used tAI.R from Mario dos Reis to compute w values for calculation of tAI
Here I will outline the steps for how I run tAI.R and how I run a permutation test to find signficiance of S-values 

The tAI.R package requires output from codonM and codonZ - I have found, however, that sometimes the order of the outputs is changed and requires "sorting" of the data 

## Step 1 - Generate the tRNA file in the order needed for tAI.R

Files required - tRNA_order.txt

This step can be done in R
```
library(tidyverse)
library(dplyr)

tRNA_file<-"INPUT_TRNA_FILE")

tRNA_data<-read_delim(tRNA_file, "\t", escape_double = FALSE, trim_ws = TRUE, col_names=TRUE))
tRNA_order<-read_delim("tRNA_order.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names=TRUE))
tRNA_new<-str_replace(tRNA_file, ".tRNA.out.tRNA.txt", ".tAI.tRNA.txt")
pec<-str_replace(tRNA_file, ".tRNA.out.tRNA.txt", "")
tRNA_data<-left_join(tRNA_order, tRNA_data)
tRNA_data<-tRNA_data[,4]
write.table(tRNA_data, file=tRNA_new, sep="", row.names = F, col.names = F)
```

You can run this first or run all the commands together 

## Step 2 - do you need to rename your sequences?

Given that the order can _sometimes_ be lost during this process I highly recommend you leave only short identifiers in your sequence names. Either the CDS or the locus tag would work in this case

for example these names contain extraneous information that might be lost therefore they should be changed to 
> >lcl|JAFVMC010000001.1_cds_KAI4938721.1_1 [locus_tag=J4E92_000001] [db_xref=InterPro:IPR010126,InterPro:IPR029058,PFAM:PF10503] [protein=hypothetical protein] [protein_id=KAI4938721.1] [location=join(7373..7717,7784..8428)] [gbkey=CDS]
> >lcl|JAFVMC010000001.1_cds_KAI4938722.1_2 [locus_tag=J4E92_000002] [db_xref=InterPro:IPR001411,InterPro:IPR011701,InterPro:IPR020846,InterPro:IPR036259,PFAM:PF06609,PFAM:PF07690] [protein=hypothetical protein] [protein_id=KAI4938722.1] [location=complement(join(8517..8662,8761..8827,9039..10316,10373..10675))] [gbkey=CDS]

## Step 3 - run codonM (see Dos Reis for more info)

```
perl codonM sequence.fna sequence.m
```

## Step 4 - run codonZ (see Dos Reis for more info)

```
codonZ sequence.fna sequence.w
```

## Step 5 - get the ws values for each codon 

This is done in R - I will load the tRNAs that I computed before


```
source("tAI.R")
library(tidyverse)
library(readr)
library(dplyr)


tRNA_file<-"FORMATTED_TRNA_FILE"
m_file<-"M_FILE.m"

#output_file_name
wi_out<-str_replace(tRNA_file, ".tAI.tRNA.txt", ".tAIR_ws.txt")

#we need the tRNA order again
suppressMessages(tRNA_order<-read_delim("tRNA_order.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names=TRUE))


this.trna <- scan(tRNA_file)
this.m <- matrix(scan(m_file), ncol=61, byrow=TRUE)
this.ws <- get.ws(tRNA=this.trna, sking=0)

#Remove these codons
ws_order<-subset(tRNA_order, AntiCodonsList != "TTA")
ws_order<-subset(ws_order, AntiCodonsList != "CTA")
ws_order<-subset(ws_order, AntiCodonsList != "TCA")
ws_order<-subset(ws_order, AntiCodonsList != "CAT")

#do some more editing - make sure the non-degenerate codons are set to 1

left_over<-subset(tRNA_order, AntiCodonsList =="TTA" |  AntiCodonsList == "CTA" |  AntiCodonsList == "TCA"|  AntiCodonsList == "CAT")
left_over$ws=1
ws_order$ws<-this.ws
ws_order<-rbind(ws_order, left_over)


ws_order<-full_join(tRNA_order,ws_order,  by="AntiCodonsList")


ws_order<-ws_order[order(ws_order$order.x),]
ws_order<-ws_order[,c(2,6)]
colnames(ws_order)<-c("CodonsList","wi")
write.table(ws_order, file=wi_out, sep="\t", row.names = F, col.names = T, quote = F)


