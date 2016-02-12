#Must use XMl package for parsing
require("XML")
#for this  we also need the foreach package for a few functions
install.packages("foreach")
library(foreach)
#for genetic distance, we use the ape package
install.packages("ape")
library(ape)

#An intial note: we want to eventually make import into R really easy by having R call the perl script directly

#original xml taken from the BOLD API based on your chosen phylum and geographical region
#xml file converted to a simplified xml with reduced number of nodes using the xmlconvertor.pl and transformxml.xsl scripts
#ensure that all files are in the same directory!
#I chose mammals as an example dataset but most likely will be switching to a smaller taxon, just doing some preliminary testing on extracting data from the nodes

#Note that taxonomic information has taxon id followed by the name of the taxon

#make sure you have the right path to the revised xml file, this path will change depending on the user!
xml <- xmlParse("C:/Users/MattL/Desktop/IGP/revisedBOLD.xml")

#Should automatically make a dataframe with the columns that we want
dfInitial <-xmlToDataFrame(getNodeSet(xml, "//record"))

#Filtering out sequences with no latitude/longitude values, filtering according to lat since we only really need lat but lon could be useful for map plotting
latlonFilter <- grep( "[0-9]", dfInitial$lat)
dfInitial<-dfInitial[latlonFilter,]

#next we have to convert lat column to num instead of chr type, this will become important later on for median latitude determination
latNum <- with(dfInitial, as.numeric(as.character(lat))) 
dfInitial$latNum <- latNum

#can do lon as well to get numeric values instead of characters for eventual map plotting
lonNum <- with(dfInitial, as.numeric(as.character(lon))) 
dfInitial$lonNum <- lonNum

#testing to see how many region nodes have numeric values since some have been reported to have coordinates in them
#region <- grep( "[0-9]", df$region)
#will tell you number of region nodes with numeric values
#length(region)

#Filtering out poor quality sequences
#As mentioned by Sally in the meeting, has a bin, N's in the middle of the sequence, greater than 500 bp and no dashes less than three in a row in the middle of a sequence = good quality sequence

#First identifying missing bins and eliminating rows with missing bin_uri's since bin is a big indicator of sequence quality
binFilter<-which(dfInitial$bin_uri=="")
if(length(binFilter) >0){
dfInitial<-dfInitial[-binFilter,]}

#To filter out sequences with N's in the middle of a sequence but not on the ends of the sequence:
nFilter <- grep( "[ACTG][N][ACTG]", dfInitial$nucleotides)
if(length(nFilter) >0){
dfInitial <- dfInitial[-nFilter,]}

#To filter out dashes less than three in the middle of a sequence since these most likely represent sequencing errors and not biologically relevant changes such as indels
dashFilter <- grep( "[ACTG][-]{1,2}[ACTG]", dfInitial$nucleotides)
if(length(dashFilter) >0){
dfInitial <- dfInitial[-dashFilter,]}

#To filter out sequences less than 500 bp
sLengths <- with(dfInitial, nchar(as.character(nucleotides))) 
dfInitial$sLengths <- sLengths
shortSeqFilter <- which(dfInitial$sLengths<500)

#finally removing the ones less than 500 bp, for this dataset
if(length(shortSeqFilter) >0){
dfInitial <- dfInitial[-shortSeqFilter,]}

#Now to separating based on latitudinal difference, starting with a default of 20 degrees latitude difference
#We first have to determine a median latitude for each bin so this will involve subdividing the dataframe into smaller dataframes by BIN

#Create groupings with each grouping representing a different bin_uri
#Each element of this list represents a bin with subelements representing the various columns of the initial dataframe created and the information is grouped by bin 
binList <- lapply(unique(dfInitial$bin_uri), function(x) dfInitial[dfInitial$bin_uri == x,])

#now to determine a median latitude for each bin
medianLat <- sapply( binList , function(x) median( x$latNum ) )

#we also need a median longitude for each if we are going to plot on a map for a visual interface
medianLon <- sapply( binList , function(x) median( x$lonNum ) )

#sapply can be used to get any statistics we want about each bin from number of members to mean of lat/lon to range of latitude in a given bin

#dataframe of our median lat values, this will be used in our final dataframe
dfLatLon <- data.frame(medianLat)

#Adding bin_uri and median longitude to dataframe
dfBin_uri <- data.frame("uri"= c(unique(dfInitial$bin_uri)))
dfLatLon$bin <- dfBin_uri
dfMedLon <- data.frame(medianLon)
dfLatLon$median <- data.frame("Lon"=c(dfMedLon$medianLon))

#Can divide into pairings based on latitude first using dfLatLon just generated, using 20 degrees latitude difference as a default
#can be found easily with dist function and put into a matrix
distLat <- dist(dfLatLon$medianLat)
distLat <- as.matrix( dist(dfLatLon$medianLat) )

#can then search matrix for values above or equal to 20 and generate a new matrix with pairings we want, give this a few seconds to load before viewing the matrix, i
distLat <- distLat[distLat[,1] >= 20,]

#we also will randomly select one sequence from each bin with all of its sequence and taxonomic data
#Foreach will iterate through each bin and sample one of its members from the initial dataframe, some of the smaller members may not be represented depending on RNG, still trying to fix this
#setting seed for testing - meaning it will produce the same result each time, if you want random each time, comment seed out
set.seed(20);
randomBinSeqList <- foreach(i=1:nrow(dfBin_uri)) %do% dfInitial[ sample( which( dfInitial$bin_uri == dfBin_uri$uri[i] ) , 1, replace = FALSE ) , ]

#Turn this list of randomly generated record ids and put into a dataframe
dfRandomBinSeq <- do.call("rbind", lapply(randomBinSeqList, as.data.frame)) 

#delete unecessary columns (region, old lat/lon etc.) of dataframe containing Random Bin members, problem is there are still a few repeated BINs so we will have to go through a few more steps to eliminate those
dfRandomBinSeq <- dfRandomBinSeq[,c(1,2,3,4,5,6,5,6,7,8,9,13)]

#Eliminates redundant bins in case there are any, now we have a usable dataframe for the distance calculation with each row being a distinct, randomly generated bin member
dfRandomBinSeq <- dfRandomBinSeq[!duplicated(dfRandomBinSeq$bin_uri),]

#Now determining sister pairs based on a genetic distance model, these pairings will then be matched with pairings generated by the latitude distance to find those that meet both criteria

#extract the nucleotide branch column from dfRandomBinSeq
dfRandomBinSeq1 <- (dfRandomBinSeq[,c(12)])

#In order to calculate genetic distances, sequences must be in a class of "DNAbin", so first we have to convert
#dfRandomBinSeq1 <-as.data.frame(dfRandomBinSeq1)
#havent figured out how to do that yet...

#Now for the computation of genetic distance, several models can be used - "raw", "N", "TS", "TV", "JC69", "K80" (the default), "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet"
#details on each model can be found here: http://www.inside-r.org/packages/cran/ape/docs/dist.dna
#Starting with a simple model called JC69 where all substituions have an equal probability
#dfGeneticDistance <- dist.dna(DNAbin, model = "JC69", as.matrix = TRUE)

#Compare function in R to find pairings that match each other:
