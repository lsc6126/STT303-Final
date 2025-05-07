# Scripts for reading in data/metadata used in analyses

###### WARNING: Lines set for use on windows,
#          Requires changes for Mac

# Gene Expression Read Counts from RNASeQCv2.4.2.
gtex_data <- read.delim(file="C:/Users/lscho/Downloads/GTEx_Analysis_2022-06-06_v10_RNASeQCv2.4.2_gene_reads.gct", skip=2, check.names=FALSE)

# Metadata for Sample Attributes
url <- "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
destfile <- "C:/Users/lscho/Downloads/GTEx_SampleAttributesDS.txt"
download.file(url, destfile, method = "curl") 
gtex_metadata <- read.delim(destfile, sep = "\t", header = TRUE)

# Metadata for Subject Phenotypes
url2 <- "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt"
destfile2 <- "C:/Users/lscho/Downloads/GTEx_SamplePhenotypesDS.txt"
download.file(url2, destfile2, method = "curl")
gtex_phenotypes <- read.delim(destfile2, sep = "\t", header = TRUE)
