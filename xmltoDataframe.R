#Must use XMl package for parsing
require("XML")

#original xml taken from the BOLD API based on your chosen phylum and geographical region
#xml file converted to a simplified xml with reduced number of nodes using the xmlconvertor.pl and transformxml.xsl scripts
#ensure that all files are in the same directory!
#I chose mammals as an example dataset

#change to your directory where the converted xml file is located
xml <- xmlParse("C:/Users/MattL/Desktop/IGP/BOLD.xml")

#Should automatically make a dataframe with the columns that we want
df <-xmlToDataFrame(getNodeSet(xml, "//record"))

#filter for rows missing species, lat/lon/. You will see afterwards that rows with missing values (except for bin_uri) have been eliminated
dfrevised<-df[complete.cases(df),]

#still a few rows with missing bin_uris
bin<-which(dfrevised$bin_uri=="")

#eliminating rows with missing bin_uri
dfrevised2<-dfrevised[-bin,]

#Filtering out poor quality sequences
#As mentioned by Sally in her email, less than 1% N (so basically no wildcards), greater than 500 bp and no dashes (i assume) = good quality sequence
#To filter out sequences with N's and dashes:
Nfilter <- grep( "[N-]", dfrevised2$nucleotides)
dfrevised3 <- dfrevised2[-Nfilter,]

#To filter out sequences less than 500 bp
sLengths <- with(dfrevised3, nchar(as.character(nucleotides))) 
dfrevised3$sLengths <- sLengths
bpshort <- which(dfrevised3$sLengths<"500")
#To check how many are less than 500 bp, if any, most seem to be around 657
length(bpshort)
#finally removing the ones less than 500 bp, for this dataset, no sequences were less than 500 bp in this dataset so we dont have to revise the dataframe again
#dfrevised4 <- dfrevised3[-bpshort,]

#Now to separating based on latitudinal difference
#Starting with a default of 20 degrees latitude difference
#We first have to determine a median latitudinal difference for each bin so this will involve subdividing the dataframe into smaller dataframes by BIN
#Create a list of dataframes with each dataframe representing a different bin_uri
bin <- lapply(unique(dfrevised3$bin_uri), function(x) dfrevised3[dfrevised3$bin_uri == x,])
#to see summary, in this case we have 92 different bins
summary(bin)
#now to determine a median latitude for each bin..
