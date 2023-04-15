#load libraries
##########################
library(tidyverse)
library(DESeq2)
library(airway)
library(readxl)
#########################

#Preparing count data
########################
setwd("/home/sam/Documents/Melanoma")
read_excel("sample.xls")