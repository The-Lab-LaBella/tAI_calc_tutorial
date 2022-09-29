# tAI_calc_tutorial
Description of my method using tAI.R (https://github.com/mariodosreis/tai)

I have used tAI.R from Mario dos Reis to compute w values for calculation of tAI
Here I will outline the steps for how I run tAI.R and how I run a permutation test to find signficiance of S-values 

The tAI.R package requires output from codonM and codonW - I have found, however, that sometimes the order of the outputs is changed and requires "sorting" of the data 

## Step 1 - Generate the tRNA file in the order needed for tAI.R

Files required - tRNA_order.txt

This step can be done in R
```
library(tidyverse)
library(dplyr)

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

Given that the order can _sometimes_ be lost during this process I highly recommend you leave only short identifiers in your sequence names

for example 
>lcl|JAFVMC010000001.1_cds_KAI4938721.1_1 [locus_tag=J4E92_000001] [db_xref=InterPro:IPR010126,InterPro:IPR029058,PFAM:PF10503] [protein=hypothetical protein] [protein_id=KAI4938721.1] [location=join(7373..7717,7784..8428)] [gbkey=CDS]
>lcl|JAFVMC010000001.1_cds_KAI4938722.1_2 [locus_tag=J4E92_000002] [db_xref=InterPro:IPR001411,InterPro:IPR011701,InterPro:IPR020846,InterPro:IPR036259,PFAM:PF06609,PFAM:PF07690] [protein=hypothetical protein] [protein_id=KAI4938722.1] [location=complement(join(8517..8662,8761..8827,9039..10316,10373..10675))] [gbkey=CDS]
>lcl|JAFVMC010000001.1_cds_KAI4938723.1_3 [locus_tag=J4E92_000003] [db_xref=InterPro:IPR017926,InterPro:IPR029062,PFAM:PF00117] [protein=hypothetical protein] [protein_id=KAI4938723.1] [location=join(10924..11152,11205..11431,11480..11770)] [gbkey=CDS]
>lcl|JAFVMC010000001.1_cds_KAI4938724.1_4 [locus_tag=J4E92_000005] [db_xref=InterPro:IPR000407,PFAM:PF01150] [protein=hypothetical protein] [protein_id=KAI4938724.1] [location=join(13093..13099,13348..15290)] [gbkey=CDS]
>
