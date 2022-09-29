# tAI_calc_tutorial
Description of my method using tAI.R (https://github.com/mariodosreis/tai)

I have used tAI.R from Mario dos Reis to compute w values for calculation of tAI
Here I will outline the steps for how I run tAI.R and how I run a permutation test to find signficiance of S-values 

The tAI.R package requires output from codonM and codonW - I have found, however, that sometimes the order of the outputs is changed and requires "sorting" of the data 

## Step 1 - Generate the tRNA file in the order needed for tAI.R

Files required - tRNA_order.txt

This step can be done in R
`library(tidyverse)
library(dplyr)

tRNA_data<-read_delim(tRNA_file, "\t", escape_double = FALSE, trim_ws = TRUE, col_names=TRUE))
tRNA_order<-read_delim("tRNA_order.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names=TRUE))
tRNA_new<-str_replace(tRNA_file, ".tRNA.out.tRNA.txt", ".tAI.tRNA.txt")
pec<-str_replace(tRNA_file, ".tRNA.out.tRNA.txt", "")
tRNA_data<-left_join(tRNA_order, tRNA_data)
tRNA_data<-tRNA_data[,4]
write.table(tRNA_data, file=tRNA_new, sep="", row.names = F, col.names = F)`
