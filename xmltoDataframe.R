#Must use XMl package for parsing
require("XML")
#For this we also need the foreach package for a few functions
install.packages("foreach")
library(foreach)
#For genetic distance, we use the ape package
install.packages("ape")
library(ape)
#Decided to use the data.tables package for one component of the code
install.packages("data.table")
library(data.table)


#An intial note: we want to eventually make import into R really easy by having R call the perl script directly
system("perl C:\\Users\\Winfield\\Desktop\\BOLD_database_xml\\xsltconverter.pl")

#Original xml taken from the BOLD API based on your chosen phylum and geographical region
#xml file converted to a simplified xml with reduced number of nodes using the xmlconvertor.pl and transformxml.xsl scripts
#Ensure that all files are in the same directory!

#Note that taxonomic information has taxon id followed by the name of the taxon
#Also note that some dataframes and matrices do not show all columns due to a limitation with R Studio 
#The command to view the full dataframe or matrix is, insert df or matrix in brackets:
#utils::View()

#Instead of pointing to the file to parse which can be tedious to constantly change the absolute path. Allows the user to pick a file and then store in a variable.
xmlParseDoc <- file.choose()
xml <- xmlParse(xmlParseDoc)

#Make sure you have the right path to the revised xml file, this path will change depending on the user!
#xml <- xmlParse("C:/Users/Winfield/Documents/revisedBOLD2.xml")

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

#Modifying Bin column slightly to remove "BIN:"
dfInitial$bin_uri <- substr(dfInitial$bin_uri, 6 , 13)

#Not filtering N's and dashes since it was mentioned by Sally that we can retain these currently

#Need to trim down each sequence if greater than 600 bp to a uniform length to be used in genetic distance calculation and in the output file
#Decided initially to trim the standard 657/658 bp sequence by 27 bp on the 5' side and 31 on the 3' side of the sequence to eliminate poor quality sequence data on either side
#Still unsure on the start position of the COI(cox) reading frame but the start and stop positions can be easily edited to match the open reading frame start of COI
dfInitial$nucleotides <- substr(dfInitial$nucleotides, 27 , 626)
dfInitial <- (dfInitial[,c("record_id","bin_uri","phylum","class","order","family","subfamily","genus","species","nucleotides","latNum","lonNum")])

#checking to make sure no sequences are not equal to 600 bp since we want a 600 bp length for the analysis
#Makes it easier as we are standardizing to one seq length intially
sLengths <- with(dfInitial, nchar(as.character(nucleotides))) 
dfInitial$sLengths <- sLengths
seqFilter <- which(dfInitial$sLength!=600)

#Finally removing the ones not equal to 600 bp
if(length(seqFilter) >0){
  dfInitial <- dfInitial[-seqFilter,]}

#We first have to determine a median latitude for each bin so this will involve subdividing the dataframe into smaller dataframes by BIN

#Create groupings by bin with each grouping representing a different bin_uri
#Each element of this list represents a bin with subelements representing the various columns of the initial dataframe created and the information is grouped by bin 
binList <- lapply(unique(dfInitial$bin_uri), function(x) dfInitial[dfInitial$bin_uri == x,])

#Now to determine a median latitude for each bin
medianLat <- sapply( binList , function(x) median( x$latNum ) )

#We also need a median longitude for each if we are going to plot on a map for a visual interface
medianLon <- sapply( binList , function(x) median( x$lonNum ) )

#we can also take a few other important pieces of data regarding each bin using sapply including number of record_ids to a bin and latitudinal min and max of each bin
latMin <- sapply( binList , function(x) min( x$latNum ) )
latMax <- sapply( binList , function(x) max( x$latNum ) )
binSize <- sapply( binList , function (x) length( x$record_id ) )

#Dataframe of our median lat values, this will be used in our final dataframe
dfLatLon <- data.frame(medianLat)

#Adding bin_uri, median longitude, latMin, latMax and binSize to dataframe with medianLat
dfLatLon$bin_uri <- c(unique(dfInitial$bin_uri))
dfLatLon$medianLon <- c(medianLon)
dfLatLon$latMin <- c(latMin)
dfLatLon$latMax <- c(latMax)
dfLatLon$binSize <- c(binSize)

#We also will randomly select one sequence from each bin of dfInitial with all of its sequence and taxonomic data, this is because for larger taxa, it may become too computionally intensive to run all sequences of each bin
#Foreach will iterate through each bin and sample one of its members from the initial dataframe, some of the smaller bins with a small number of members may not necessarily be represented depending on RNG (random number generation) but if run enough times should represent all BINs in the chosen taxa
#Had set a seed for testing - meaning it will produce the same result each time, if you want random each time, comment seed out
set.seed(15);
randomBinSeqList <- foreach(i=1:nrow(dfLatLon)) %do% dfInitial[ sample( which( dfInitial$bin_uri == dfLatLon$bin_uri[i] ) , 1, replace = FALSE ) , ]

#Turn this list of randomly sampled record ids and put into a dataframe
dfRandomBinSeq <- do.call("rbind", lapply(randomBinSeqList, as.data.frame)) 

#Eliminates redundant bins in case there are any, now we have a usable dataframe for the distance calculation with each row being a distinct, randomly generated bin member
dfRandomBinSeq <- dfRandomBinSeq[!duplicated(dfRandomBinSeq$bin_uri),]

#Can reference to our median Latitude dataframe and determine which ones are shared since some bins might not necessarily be represented due to the sample function used earlier
dfLatLon <- subset(dfLatLon, dfLatLon$bin_uri %in% dfRandomBinSeq$bin_uri)

#Can now merge dfLatLon with RandomBinseq (but retain the name RandomBinSeq) to create one dataframe containg all pertinent data we need and can reference back to from the match overall dataframe produced later on
dfRandomBinSeq <- merge(dfRandomBinSeq, dfLatLon, by.x = "bin_uri")
#Getting rid of latNum and lonNum since we dont need these anymore, we only need the median Lat/Lon values 
dfRandomBinSeq <- (dfRandomBinSeq[,c("bin_uri","binSize","record_id","phylum","class","order","family","subfamily","genus","species","nucleotides","sLengths","medianLat","latMin","latMax","medianLon")])
#Adding an index column to reference later with Match overall dataframe
dfRandomBinSeq$ind <- row.names(dfRandomBinSeq)

#Can then determine latitude differences, originally put into a matrix but transformed into a df
#Can be found easily with dist function and put into a matrix
distLat <- dist(dfRandomBinSeq$medianLat)
matrixLatitudeDistance <- as.matrix( dist(dfRandomBinSeq$medianLat) )
#Convert to dataframe since its easier to manipulate, this dataframe will be used further down
dfLatitudeDistance <-as.data.frame(matrixLatitudeDistance)

#Now determining sister pairs based on a genetic distance model, these pairings will then be matched with pairings generated by the latitude distance to find those that meet both criteria
#Extracting the nucleotide column from dfRandomBinSeq to grab the actual sequencing data
dfRandomBinSeqCol <- (dfRandomBinSeq[,c("nucleotides")])

#In order to calculate genetic distances, sequences must be in a class of "DNAbin" according to the Ape package, so first we have to split by base pair, align the sequences and then convert to DNAbin format for use in the Ape package
randomBinSeqSplit <-strsplit(dfRandomBinSeqCol, "")
DNAalign <- as.alignment(randomBinSeqSplit)
DNAbin <- as.DNAbin(DNAalign)

#Now for the computation of genetic distance, several models can be used - "raw", "N", "TS", "TV", "JC69", "K80" (the default), "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet"
#Details on each model can be found here: http://www.inside-r.org/packages/cran/ape/docs/dist.dna
#Starting with the K80 model since this was one of the models suggested by Sally
#Note that for this comparison all sequences must be equal in length but this is already taken care of at an earlier step
matrixGeneticDistance <- dist.dna(DNAbin, model = "K80", as.matrix = TRUE, pairwise.deletion = TRUE)
#convert to dataframe)
dfGeneticDistance <-as.data.frame(matrixGeneticDistance)
#Putting it into a stack (each column concatenated into one long column of indexes and values) so it can be easily subsetted
dfGeneticDistanceStack <-stack(dfGeneticDistance)

#So now we have two dataframes, we can find ideal matchings based on less than 15% divergence and 20 degrees latitude separation
#These values can easily be edited to suit your needs
#Will produce lists with indexes of each match according to our set criteria
geneticDistanceMatchI <- which(dfGeneticDistance<=0.15)
latitudeDistanceMatch <- which(dfLatitudeDistance>=20)

#Match the dfs against each other to find the sister pairs using intersect
matchOverall <- intersect(geneticDistanceMatchI,latitudeDistanceMatch)

#Use these index matches and reference against the genetic distance stack dataframe to subset it according to matches
dfMatchOverall <- dfGeneticDistanceStack[c(matchOverall), ]
#Each duplication of genetic distance values corresponds to a pairing 
#Will now merge the randombinseq dataframe to the pairings generated 
dfMatchOverall <- merge(dfRandomBinSeq, dfMatchOverall, by.x = "ind", by.y = "ind")

#Pairings will be ordered according to genetic distance, least divergent pairing to most divergent pairing
dfMatchOverall <- dfMatchOverall[order(dfMatchOverall$values),] 

#Then we can multiply these distance values by 1.3 to determine the ideal minimum outgroup distance from the pairings, we put this in another column in the matchOverall dataframe
#This minimum outgroup distance could also be a user adjustable parameter
dfMatchOverall$inGroupDistx1.3 <- dfMatchOverall$values * 1.3

#Also grouping dataframe every 2 rows to reflect each unique pairing/matching, pairing column will give a number value for each pairing ordered 
dfMatchOverall$inGroupPairing <- rep(1:(nrow(dfMatchOverall)/2), each = 2)

#Reorganizing and renaming some columns in MatchOverall dataframe to make more easily readable 
dfMatchOverall <- (dfMatchOverall[,c("inGroupPairing","bin_uri","record_id","values","inGroupDistx1.3","medianLat","latMin","latMax","binSize","phylum","class","order","family","subfamily","genus","species","nucleotides","sLengths","medianLon","ind")])
colnames(dfMatchOverall)[4] <- "inGroupDist"
colnames(dfMatchOverall)[6] <- "binMedianLat"
colnames(dfMatchOverall)[7] <- "binLatMin"
colnames(dfMatchOverall)[8] <- "binLatMax"
colnames(dfMatchOverall)[10] <- "taxon_id:phylum"
colnames(dfMatchOverall)[11] <- "taxon_id:class"
colnames(dfMatchOverall)[12] <- "taxon_id:order"
colnames(dfMatchOverall)[13] <- "taxon_id:family"
colnames(dfMatchOverall)[14] <- "taxon_id:subFamily"
colnames(dfMatchOverall)[15] <- "taxon_id:genus"
colnames(dfMatchOverall)[16] <- "taxon_id:species"
colnames(dfMatchOverall)[17] <- "sequenceCOI"
colnames(dfMatchOverall)[20] <- "indexNo"

#For bins that have multiple pairings we can determine which pairing has the smallest divergence and remove the others in order to choose the best possible one
#Basically we are taking the best possible pairing per bin based on divergence
#I will call this one dfMatchOverall1
dfMatchOverall1 <- by(dfMatchOverall, dfMatchOverall["bin_uri"], head, n=1)
dfMatchOverall1 <- Reduce(rbind, dfMatchOverall1)
dfMatchOverall1 <- subset(dfMatchOverall1, duplicated(dfMatchOverall1$inGroupPairing))
#Can then take our best pairings and subset them against the original dfMatchOverall dataframe, I will call this dfMatchOverallBest
dfMatchOverallBest <- subset(dfMatchOverall, dfMatchOverall$inGroupPairing %in% dfMatchOverall1$inGroupPairing)

#We now have the best possible pairings from each bin based on both latitude and distance criteria!

#We now have to determine the best possible outgroupings

#Can find suitable outgroupings based on the calculated mimumum genetic distances of the ingroupdist1.3x column and indexNo columns
#First we can use the indexNo to subset the dfGeneticDistanceStackDataframe according to indexes represented in the final pairings, this will be called dfBestOutGroup
#This will essentially limit our outgroup bin distances to those associated with the ingroup bin distances
dfBestOutGroup <- subset(dfGeneticDistanceStack, dfGeneticDistanceStack$ind %in% dfMatchOverall1$indexNo)

#Determining shared indices and putting in a variable
outGroupCandidates <- foreach(i=1:nrow(dfMatchOverall1)) %do% which(dfBestOutGroup$ind == dfMatchOverall1$indexNo[i])
#Determining indices with correct outgroup distance 
outGroupCandidates1 <- foreach(i=1:nrow(dfMatchOverall1)) %do% which(dfBestOutGroup$values >= dfMatchOverall1$inGroupDistx1.3[i])
#Intersection of the two using mapply to find the correct outgroupings for each pairing
outGroupCandidates2 <- mapply(intersect,outGroupCandidates,outGroupCandidates1)
#Unlist to make one easy to handle vector
outGroupCandidates3 <- unlist(outGroupCandidates2)
#Adding an nrow column to dfBestOutGroup for subsetting
dfBestOutGroup$rownum<-seq.int(nrow(dfBestOutGroup))
#Then we can subset to nrow column based on the outgroup candidates
dfBestOutGroup <- dfBestOutGroup[dfBestOutGroup$rownum %in% outGroupCandidates3, ]

#Break this down into a list of elements, each element being a unique index
outGroupList <- lapply(unique(dfBestOutGroup$ind), function(x) dfBestOutGroup[dfBestOutGroup$ind == x,])
#Use sapply to find the min value of each list element, the min value being the closest value to the outgroup dist we want (x1.3ingroupdist)
closestOutGroup <- sapply( outGroupList , function(x) min( x$values ) )
#Turn into data.frame
dfClosestOutGroup <- as.data.frame(closestOutGroup)
#Again subset to dfBestOutGroup to get indices of the closest outgrouping values
dfBestOutGroup <- dfBestOutGroup[dfBestOutGroup$values %in% dfClosestOutGroup$closestOutGroup,]

#Its possible to have more than one ideal outgrouping per pairing, so restricting to one outgrouping per pairing
dfBestOutGroup <- by(dfBestOutGroup, dfBestOutGroup["ind"], head, n=1)
dfBestOutGroup <- Reduce(rbind, dfBestOutGroup)
dfBestOutGroup<- dfBestOutGroup[order(dfBestOutGroup$rownum),] 

#More Dataframe/table manipulation

#We also need the identity of our outgroupings, this will give numerical identity in order to reference dfrandombinseq to get outgrouping species data
j=0:(nrow(dfBestOutGroup)-1)
dfBestOutGroup$identity <- (dfBestOutGroup$rownum / (nrow(dfRandomBinSeq))-j) * (nrow(dfRandomBinSeq))
#Also ordering dfBestOutGroup by identity
dfBestOutGroup<- dfBestOutGroup[order(dfBestOutGroup$identity),] 
colnames(dfRandomBinSeq)[17] <- "identity"
#Converting to data table format and setting keys in order to properly merge
dfRandomBinSeq <- data.table(dfRandomBinSeq)
dfBestOutGroup <- data.table(dfBestOutGroup)
setkey(dfBestOutGroup,identity)
#convert identity key of dfRandomBinSeq to numeric
identityNum <- with(dfRandomBinSeq, as.numeric(as.character(identity))) 
dfRandomBinSeq$identityNum <- identityNum
#Reorder dfRandomBinSeq and setting key for identityNum
dfRandomBinSeq<- dfRandomBinSeq[order(dfRandomBinSeq$identityNum),] 
setkey(dfRandomBinSeq,identityNum)

#Can now merge randombinseq to bestoutgroup giving us the data we need
dfBestOutGroup <- merge(dfBestOutGroup, dfRandomBinSeq, by.x = "identity", by.y = "identityNum", all.x = TRUE) 
dfBestOutGroup <- as.data.frame(dfBestOutGroup)

#Can rename certain columns to more closely resemble dfMatchOverall 
dfBestOutGroup <- dfBestOutGroup[,c("bin_uri","record_id","values","medianLat","latMin","latMax","binSize","phylum","class","order","family","subfamily","genus","species","nucleotides","sLengths","medianLon","identity")]
colnames(dfBestOutGroup)[3] <- "outGroupDist"
colnames(dfBestOutGroup)[4] <- "binMedianLat"
colnames(dfBestOutGroup)[5] <- "binLatMin"
colnames(dfBestOutGroup)[6] <- "binLatMax"
colnames(dfBestOutGroup)[8] <- "taxon_id:phylum"
colnames(dfBestOutGroup)[9] <- "taxon_id:class"
colnames(dfBestOutGroup)[10] <- "taxon_id:order"
colnames(dfBestOutGroup)[11] <- "taxon_id:family"
colnames(dfBestOutGroup)[12] <- "taxon_id:subFamily"
colnames(dfBestOutGroup)[13] <- "taxon_id:genus"
colnames(dfBestOutGroup)[14] <- "taxon_id:species"
colnames(dfBestOutGroup)[15] <- "sequenceCOI"
colnames(dfBestOutGroup)[18] <- "indexNo"

#Now just need to associate our outgroup with the pairings by adding another column to the outgroup df 
#Just ordering values by order of values in ClosestOutGroup to ensure we have the correct order
dfBestOutGroup <- dfBestOutGroup[order(match(dfBestOutGroup[,3],dfClosestOutGroup[,1])),]

#Appending the ingrouppairing column to the Outgroup dataframe so we know which pairing its associated with
dfBestOutGroup$associatedInGroup <- dfMatchOverall1$inGroupPairing
#One more reorganization
dfBestOutGroup <- (dfBestOutGroup[,c("associatedInGroup","bin_uri","record_id","outGroupDist","binMedianLat","binLatMin","binLatMax","binSize","taxon_id:phylum","taxon_id:class","taxon_id:order","taxon_id:family","taxon_id:subFamily","taxon_id:genus","taxon_id:species","sequenceCOI","sLengths","medianLon","indexNo")])

#Can now merge eveything together into the MatchOverallBest, this will act as our finalized dataframe
#Note that each outgroup will be duplicated for each member of the pairing
#.x beside headings are the ingroup data (except unique columns), scroll far enough to the right and headings with .y are the outgroup related columns
dfMatchOverallBest <- merge(dfMatchOverallBest, dfBestOutGroup, by.x = "inGroupPairing", by.y = "associatedInGroup")

#Note there will be a few pairings that dont have a suitable outgroup(couldnt find a bin with the suitable outgroup distance) so we can filter those out
colnames(dfMatchOverallBest)[21] <- "bin_uriOutGroup"
binFilter2<-which(is.na(dfMatchOverallBest$bin_uriOutGroup))
if(length(binFilter2) >0){
dfMatchOverallBest<-dfMatchOverallBest[-binFilter2,]}

#Now we have both outgroups and suitable pairings!

#Ouput the contents to a CSV file in CSV format

