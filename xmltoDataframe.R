#Must use XMl package for parsing
require("XML")
#For this we also need the foreach package for a few functions
install.packages("foreach")
library(foreach)
#For genetic distance, we use the ape package
install.packages("ape")
library(ape)

#An intial note: we want to eventually make import into R really easy by having R call the perl script directly

#Original xml taken from the BOLD API based on your chosen phylum and geographical region
#xml file converted to a simplified xml with reduced number of nodes using the xmlconvertor.pl and transformxml.xsl scripts
#Ensure that all files are in the same directory!
#I chose mammals as an example dataset but most likely will be switching to a smaller taxon, just doing some preliminary testing on extracting data from the nodes

#Note that taxonomic information has taxon id followed by the name of the taxon

#Make sure you have the right path to the revised xml file, this path will change depending on the user!
xml <- xmlParse("C:/Users/MattL/Desktop/IGP/revisedBOLD.xml")

#Should automatically make a dataframe with the columns that we want
dfInitial <-xmlToDataFrame(getNodeSet(xml, "//record"))

#Removing sequences with no latitude/longitude values, filtering according to lat since we only really need lat but lon could be useful for map plotting
containLatLon <- grep( "[0-9]", dfInitial$lat)
dfInitial<-dfInitial[containLatLon,]

#Next we have to convert lat column to num instead of chr type, this will become important later on for median latitude determination
latNum <- with(dfInitial, as.numeric(as.character(lat))) 
dfInitial$latNum <- latNum

#Can do lon as well to get numeric values instead of characters for eventual map plotting
lonNum <- with(dfInitial, as.numeric(as.character(lon))) 
dfInitial$lonNum <- lonNum

#Filtering out poor quality sequences
#First identifying missing bins and eliminating rows with missing bin_uri's since bin is a big indicator of sequence quality
binFilter<-which(dfInitial$bin_uri=="")
if(length(binFilter) >0){
dfInitial<-dfInitial[-binFilter,]}

#Not filtering N's and dashes since it was mentioned that we can retain these currently, most of these will be removed in the trimming anyways

#Need to trim down each sequence if greater than 600 bp to a uniform length to be used in genetic distance calculation and in the output file
#Decided initially to trim the standard 657/658 bp sequence by 27 bp on the 5' side and 31 on the 3' side of the sequence to eliminate poor quality sequence data on either side
#Still unsure on the start position of the COI(cox) reading frame but the start and stop positions can be easily edited to match the open reading frame start of COI
Trim <- substr(dfInitial$nucleotides, 27 , 626)
dfTrim <- data.frame(Trim)
#Append this column to dfInitial and delete old nucleotides and lat/lon columns
dfInitial$seq <- Trim
dfInitial <- (dfInitial[,c("record_id","bin_uri","phylum","class","order","family","subfamily","genus","species","seq","latNum","lonNum")])

#Checking to make sure no sequences are not equal to 600 bp since we want a 600 bp length for the analysis
#Makes it easier as we are standardizing to one seq length intially
sLengths <- with(dfInitial, nchar(as.character(seq))) 
dfInitial$sLengths <- sLengths
seqFilter <- which(dfInitial$sLength!=600)

#Finally removing the ones not equal to 600 bp
if(length(seqFilter) >0){
  dfInitial <- dfInitial[-seqFilter,]}

#Create groupings by bin with each grouping representing a different bin_uri
#Each element of this list represents a bin with subelements representing the various columns of the initial dataframe created and the information is grouped by bin 
binList <- lapply(unique(dfInitial$bin_uri), function(x) dfInitial[dfInitial$bin_uri == x,])

#Now to determine a median latitude for each bin
medianLat <- sapply( binList , function(x) median( x$latNum ) )

#We also need a median longitude for each if we are going to plot on a map for a visual interface
medianLon <- sapply( binList , function(x) median( x$lonNum ) )

#sapply can be used to get any statistics we want about each bin from number of members to mean of lat/lon to range of latitude in a given bin 

#Dataframe of our median lat values, this will be used in our final dataframe
dfLatLon <- data.frame(medianLat)

#Adding bin_uri and median longitude to dataframe
dfBin_uri <- data.frame("uri"= c(unique(dfInitial$bin_uri)))
dfLatLon$bin <- dfBin_uri
dfMedLon <- data.frame(medianLon)
dfLatLon$median <- data.frame("Lon"=c(dfMedLon$medianLon))

#We also will randomly select one sequence from each bin of dfInitial with all of its sequence and taxonomic data
#Foreach will iterate through each bin and sample one of its members from the initial dataframe, some of the smaller members may not be represented depending on RNG, still trying to fix this
#Had set a seed for testing - meaning it will produce the same result each time, if you want random each time, comment seed out
set.seed(10);
randomBinSeqList <- foreach(i=1:nrow(dfBin_uri)) %do% dfInitial[ sample( which( dfInitial$bin_uri == dfBin_uri$uri[i] ) , 1, replace = FALSE ) , ]

#Turn this list of randomly sampled record ids and put into a dataframe
dfRandomBinSeq <- do.call("rbind", lapply(randomBinSeqList, as.data.frame)) 

#Eliminates redundant bins in case there are any, now we have a usable dataframe for the distance calculation with each row being a distinct, randomly generated bin member
dfRandomBinSeq <- dfRandomBinSeq[!duplicated(dfRandomBinSeq$bin_uri),]

#Can reference to our median Latitude dataframe and determine which ones are shared
dfLatLon <- subset(dfLatLon, dfBin_uri$uri %in% dfRandomBinSeq$bin_uri)

#Can then determine latitude differences, originally put into a matrix but transformed into a df
#Can be found easily with dist function and put into a matrix
distLat <- dist(dfLatLon$medianLat)
matrixLatitudeDistance <- as.matrix( dist(dfLatLon$medianLat) )
#Convert to dataframe since its easier to manipulate, this dataframe will be used further down
dfLatitudeDistance <-as.data.frame(matrixLatitudeDistance)

#Now determining sister pairs based on a genetic distance model, these pairings will then be matched with pairings generated by the latitude distance to find those that meet both criteria
#Extract the nucleotide column from dfRandomBinSeq
dfRandomBinSeqCol <- (dfRandomBinSeq[,c("seq")])

#In order to calculate genetic distances, sequences must be in a class of "DNAbin", so first we have to split by base pair, align the sequences and then convert to DNAbin format for use in the Ape package
randomBinSeqSplit <-strsplit(dfRandomBinSeqCol, "")
DNAalign <- as.alignment(randomBinSeqSplit)
DNAbin <- as.DNAbin(DNAalign)

#Now for the computation of genetic distance, several models can be used - "raw", "N", "TS", "TV", "JC69", "K80" (the default), "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet"
#Details on each model can be found here: http://www.inside-r.org/packages/cran/ape/docs/dist.dna
#Starting with the K80 model since this was one of the models suggested by Sally
#Note that for this comparison all sequences must be equal in length but this is already taken care of at an earlier step
matrixGeneticDistance <- dist.dna(DNAbin, model = "K80", as.matrix = TRUE)
#convert to dataframe
dfGeneticDistance <-as.data.frame(matrixGeneticDistance)
#Putting it into a stack (each column concatenated into one long column of indexes and values) so it can be easily subsetted
dfGeneticDistanceStack <-stack(dfGeneticDistance)

#So now we have two dataframes - one with lat distance and the other with genetic distance, we can find ideal matchings based on less than 15% divergence and 20 degrees latitude separation
#These values can easily be edited to suit your needs
#Will produce dataframes with indexes of each match according to our set criteria
geneticDistanceMatch <- which(dfGeneticDistance<=0.15)
latitudeDistanceMatch <- which(dfLatitudeDistance>=20)

#Match the dfs against each other to find the sister pairs using intersect
matchOverall <- intersect(geneticDistanceMatch,latitudeDistanceMatch)
#Use these index matches and reference against the genetic distance dataframe to subset it according to matches
dfMatchOverall <- dfGeneticDistanceStack[c(matchOverall), ]

#Then we can multiply these distance values by 1.3 to determine the minimum outgroup distance from the pairings, we put this in another column in the matchOverall dataframe
dfMatchOverall$valuesx1.3 <- dfMatchOverall$values * 1.3


#Now have to determine outgroupings...



#Once matching pairs plus outgroups have been determined, can attach all relevant data to those in dataframe




#Ouput the dataframe to a CSV file in CSV format
