#################
#Authored by BIF IGP Group Tempura
#Major contributions by Matthew Orton to this script

#This program will allow for the generation of latitudinally separated sister pairings and associated outgroupings from ANY taxa and geographical region found on BOLD in one streamlined R pipeline!
#In this new and improved iteration, data is translated directly from a BOLD tsv file of the users choosing to a dataframe in R
#The generated sister pairs and outgroups can then be written to a csv or tsv and the file will appear in the current working directory of R

##################
#A few important tips:

#There are two options for parsing the tsv, you can either download the tsv directly from the BOLD API OR you can use a tsv you have previously downloaded

#Larger taxa, ex: aves can potentially take from several minutes to hours in order to process the tsv and run the sequence alignment. This can potentially consume a lot of working memory. I wouldnt run any taxa that is very large, (insecta for example) until we know what kind of memory resources it will consume.

#Some tips for using the BOLD API since this is what is used to grab the relevant data we need: 
#To see details on how to use the bold API, go to http://www.boldsystems.org/index.php/resources/api?type=webservices
#Can add additional restrictions to url, for example &instituiton=Biodiversity Institute of Ontario|York University or &marker=COI-5P if you want to specifiy an institution or specific genetic marker in the url
#geo=all in the url means global but geo can be geo=Canada for example
#Can use | modifier in the url, for example &geo=Canada|Alaska would give data for both Canada and alaska, or taxon=Aves|Reptilia would yield give data for both Aves and Reptilia

#We recommend using RStudio: https://www.rstudio.com/home/, provides a nicer interface than just using R
#Large matrices and dataframes may not show all columns, you can overcome this by typing the command: utils::View() where you would insert the dataframe or matrix you want in the brackets of that command

#Important dataframes:
#dfMatchOverallBest is the finalized dataframe that contains all of the finalized pairings and outgroupings
#dfMatchOverall is the dataframe of ingroup pairings before bins with multiple pairings are eliminated based on distance and before the problem of pseudoreplication is addressed
#dfMatchOverallLineage1 and Lineage2 represent dataframes each member of the pairing
#dfInitial is the dataframe first produced by the import from BOLD and is trimmed by lat, bin_uri etc.
#dfLatLon just contains relevant information for each bin: bin size, maximum lat, median lat, minimum lat, median lon
#dfBestOutGroupL1 and L2 contain the associated outgroupings only (for each lineage) but each does have a column for which pairing its associated with
#dfLatitudeDistance and dfGeneticDistance are matrices converted into dataframes that show all possible distances between bins
#dfAllSeq is the dataframe that contains sequence data for both consensus sequences and nonconsensus (bins with one member to them) sequences
#dfGeneticDistanceStack is dfGeneticDistance with all columns concatenated into one long column, it is used to grab index numbers for each pairing
#dfDistancePair represents pairwise distances between lineages in the final pairings only

#################
#Packages required

#For this we also need the foreach package for a few functions
install.packages("foreach")
library(foreach)
#For genetic distance, we use the ape package
install.packages("ape")
library(ape)
#Speeds up parsing of the tsv file with read_tsv function
install.packages("readr")
library(readr)
#For sequence alignments we need the biostrings (DNAStringSet function) and msa packages, run each of these commands individually, sometimes it skips the libraries
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
library("Biostrings")
biocLite("msa")
library("msa")
#################
#R Commands:

#TSV Parsing

#First we download the TSV and convert it into a dataframe, this URL is what is modified by the user and will determine the taxa, geographic region etc.
dfInitial <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Tardigrada&geo=all&format=tsv")

#If you want to run pre downloaded BOLD tsv's, this will let you choose a path to that tsv and parse
#tsvParseDoc <- file.choose()
#dfInitial <- read_tsv(tsvParseDoc)

##############
#Dataframe Filtering and Reorganization

#Filtering this df according to the relevant columns we need
dfInitial <- (dfInitial[,c("recordID","bin_uri","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","lat","lon","nucleotides")])
colnames(dfInitial)[1] <- "record_id"

#Removing sequences with no latitude/longitude values, filtering according to lat since we only really need lat but lon could be useful for map plotting
containLatLon <- grep( "[0-9]", dfInitial$lat)
dfInitial<-dfInitial[containLatLon,]

#Next we have to convert lat column to num instead of chr type, this will become important later on for median latitude determination
latNum <- with(dfInitial, as.numeric(as.character(lat))) 
dfInitial$latNum <- latNum

#Can do lon as well to get numeric values instead of characters
lonNum <- with(dfInitial, as.numeric(as.character(lon))) 
dfInitial$lonNum <- lonNum

#First identifying missing bins and eliminating rows with missing bin_uri's since bin is a big indicator of sequence quality
#Grep by colon since every record with a bin identifier will have this
containBin <- grep( "[:]", dfInitial$bin_uri)
dfInitial<-dfInitial[containBin,]

#Getting rid of any records that dont have sequence data (sometimes there are a few)
containNucleotides <- grep( "[ACGTN]", dfInitial$nucleotides)
dfInitial<-dfInitial[containNucleotides,]

#Modifying Bin column slightly to remove "BIN:"
dfInitial$bin_uri <- substr(dfInitial$bin_uri, 6 , 13)

#Dataframe Reorganization
dfInitial <- (dfInitial[,c("record_id","bin_uri","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","latNum","lonNum")])

############
#Bin Stats and Median Latitude/Longitude Determination per bin

#First we can make a smaller dataframe with the columns we want for each bin - bin_uri, latnum, lonnum, record_id, if we didnt do this, the binList would consume a huge amount of memory
dfBinList <- (dfInitial[,c("record_id","bin_uri","latNum","lonNum","nucleotides")])
#Create groupings by bin with each grouping representing a different bin_uri
#Each element of this list represents a bin with subelements representing the various columns of the initial dataframe created and the information is grouped by bin 
binList <- lapply(unique(dfBinList$bin_uri), function(x) dfBinList[dfBinList$bin_uri == x,])

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

#Can also convert to 180 degree latitude scale, convert longitude as well
dfLatLon$medianLat <- dfLatLon$medianLat + 90
dfLatLon$latMax <- dfLatLon$latMax + 90
dfLatLon$latMin <- dfLatLon$latMin + 90

#Merging LatLon to BinList for the sequence alignment step
dfBinList <- merge(dfBinList, dfLatLon, by.x = "bin_uri", by.y = "bin_uri")

###############
#Sequence alignments and consensus sequences for each bin with bins larger than 1 member

#First we have to subset dfBinList to find bins with more than one member since these bins will need consensus sequences
largeBin <-which(dfBinList$binSize > 1)
#If there is at least one bin with more than one member, then a dataframe dfConsensus will be created with those bins
if(length(largeBin) >0){
  dfConsensus <- dfBinList[largeBin,]

  #Also need to find the number of unique bins in dfConsensus
  binNumberConsensus <- unique(dfConsensus$bin_uri)
  binNumberConsensus <- length(binNumberConsensus)

  #We also have to create another separate dataframe with bins that only have one member
  dfNonConsensus <- dfBinList[-largeBin,]

  #We then take the dfConsensus sequences and break it down into a list with each element being a unique bin
  largeBinList <- lapply(unique(dfConsensus$bin_uri), function(x) dfConsensus[dfConsensus$bin_uri == x,])
  
  #Convert all of the sequences in this list to dnaStringSet format for the alignment step
  dnaStringSet1 <- sapply( largeBinList, function(x) DNAStringSet(x$nucleotides) )

  #Multiple sequence alignment using msa package - 3 different algorithms can be used for this: ClustalW, ClustalOmega, MUSCLE, using the default of ClustalW for testing
  #Sequence alignment package can be found here: http://master.bioconductor.org/packages/3.2/bioc/vignettes/msa/inst/doc/msa.pdf
  #Using the default parameters of the msa package (ClustalW)
  #Run a multiple sequence alignment on each element of the dnaStringSet1 list, can achieve this with the foreach package
  alignment1 <- foreach(i=1:binNumberConsensus) %do% msa(dnaStringSet1[[i]])
  
  #We can then extract the consensus sequence from each individual alignment
  consensusSeq <- foreach(i=1:binNumberConsensus) %do% consensusString(alignment1[[i]])
  
  #Adding out consensusSeq to dfConsensus 
  dfConsensus <- by(dfConsensus, dfConsensus["bin_uri"], head, n=1)
  dfConsensus <- Reduce(rbind, dfConsensus)
  dfConsensus$nucleotides <- consensusSeq
  #Labeling all consensus sequences as such in the record_id since they would not have a specific record_id
  dfConsensus$record_id <- "consensus sequence"
  
  #Append our nonconsensus to our consensus, will make a new dataframe for this - dfAllSeq
  #now we have all of the sequences we need for the next alignment of all sequences
  dfAllSeq <-  rbind(dfConsensus, dfNonConsensus)
  #Can then merge this with dfInitial to get all of the relevant data we need
  #merging to dfInitial to gain all the relevant taxanomic data
  #of course consensus sequences will not have any specific taxonomic data but will still retain their bin_uri for identification
  dfAllSeq <- merge(dfAllSeq, dfInitial, all.x = TRUE)
  #Renaming and reorganizing the dataframe
  dfAllSeq <- (dfAllSeq[,c("bin_uri","binSize","record_id","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","medianLat","latMin","latMax","medianLon")])
  #Adding an index column to reference later with Match overall dataframe
  dfAllSeq$ind <- row.names(dfAllSeq)
  
  } else {
  #Else if there are no bins with more than one member than we would simply merge latlon with intial to get dfAllSeq and thus all of the sequences we want
  dfAllSeq <- merge(dfLatLon, dfInitial, by.x = "record_id", by.y = "record_id")
  dfAllSeq <- (dfAllSeq[,c("bin_uri","binSize","record_id","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","medianLat","latMin","latMax","medianLon")])
  dfAllSeq$ind <- row.names(dfAllSeq)
}

#############
#Latitude Distance Determination

#Can be determined easily with the dist function which will determine latitudinal differences between all bins
distLat <- dist(dfAllSeq$medianLat)
#Then we convert this to a matrix
matrixLatitudeDistance <- as.matrix( dist(dfAllSeq$medianLat) )
#Then convert to dataframe since its easier to manipulate, this dataframe will be used further down
dfLatitudeDistance <-as.data.frame(matrixLatitudeDistance)

#############
#Multiple Sequence Alignment of All Sequences and Pairwise Distance determination with TN93
#This involves aligning all sequences and using those sequences in a pairwise distance computation

#First have to convert the nucleotide column to type character for this to work
dfAllSeq$nucleotides <- with(dfAllSeq, as.character(nucleotides)) 

#Converting all sequences to DNAStringSet format again
dnaStringSet2 <- DNAStringSet(dfAllSeq$nucleotides)

#Run a multiple sequence alignment of all sequences both consensus and nonconsensus
#As mentioned using the default of ClustalW
#This could take several minutes to even hours in some cases depending on the taxa
alignment2 <- msa(dnaStringSet2)

#Conversion to DNAbin format before genetic distance computation matrix
DNAbin <- as.DNAbin(alignment2)

#Now for the computation of genetic distance, several models can be used - "raw", "N", "TS", "TV", "JC69", "K80" (the default), "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet"
#Details on each model can be found here: http://www.inside-r.org/packages/cran/ape/docs/dist.dna
#Using the TN93 model for our data
matrixGeneticDistance <- dist.dna(DNAbin, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)
#convert to dataframe
dfGeneticDistance <-as.data.frame(matrixGeneticDistance)
#Putting it into a stack (each column concatenated into one long column of indexes and values) so it can be easily subsetted
dfGeneticDistanceStack <-stack(dfGeneticDistance)

##############
#Finding appropriate pairings according to latitude and distance criteria

#So now we have two dataframes, we can find ideal matchings based on less than 15% divergence and 20 degrees latitude separation
#These values can easily be edited to suit the user
#Will produce lists with indexes of each match according to our set criteria
geneticDistanceMatchI <- which(dfGeneticDistance<=0.15)
latitudeDistanceMatch <- which(dfLatitudeDistance>=20)

#Match the dfs against each other to find the sister pairs using intersect
matchOverall <- intersect(geneticDistanceMatchI,latitudeDistanceMatch)

#Use these index matches and reference against the genetic distance stack dataframe to subset it according to matches
dfMatchOverall <- dfGeneticDistanceStack[c(matchOverall), ]
#Each duplication of genetic distance values corresponds to a pairing 
#Will now merge the randombinseq dataframe to the pairings generated 
dfMatchOverall <- merge(dfAllSeq, dfMatchOverall, by.x = "ind", by.y = "ind")

#Pairings will be ordered according to genetic distance, least divergent pairing to most divergent pairing
dfMatchOverall <- dfMatchOverall[order(dfMatchOverall$values),] 

#Then we can multiply these distance values by 1.3 to determine the ideal minimum outgroup distance from the pairings, we put this in another column in the matchOverall dataframe
#This minimum outgroup distance could also be a user adjustable parameter
dfMatchOverall$inGroupDistx1.3 <- dfMatchOverall$values * 1.3

#Also grouping dataframe every 2 rows to reflect each unique pairing/matching, pairing column will give a number value for each pairing ordered 
dfMatchOverall$inGroupPairing <- rep(1:(nrow(dfMatchOverall)/2), each = 2)

#Reorganizing and renaming some columns in MatchOverall dataframe to make more easily readable 
dfMatchOverall <- (dfMatchOverall[,c("inGroupPairing","record_id","bin_uri","values","inGroupDistx1.3","medianLat","latMin","latMax","binSize","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","medianLon","ind")])
colnames(dfMatchOverall)[4] <- "inGroupDist"
colnames(dfMatchOverall)[26] <- "indexNo"

##############
#Taking the Best Possible Pairings and Creating Dataframes for each Lineage and each Complete Pairing

#For bins that have multiple pairings associated with them we can select which pairing has the smallest divergence to its other ingroup member since the pairings are ordered according to distance
#Basically we are taking the best possible pairing per bin based on divergence such that a bin would never be duplicated in multiple pairings
#This dataframe will be dfMatchOverallLineage1
dfMatchOverallLineage1 <- by(dfMatchOverall, dfMatchOverall["bin_uri"], head, n=1)
#This dataframe represents the first lineage of each sister pairing
dfMatchOverallLineage1 <- Reduce(rbind, dfMatchOverallLineage1)
dfMatchOverallLineage1 <- subset(dfMatchOverallLineage1, duplicated(dfMatchOverallLineage1$inGroupPairing))

#To get both lineages for each pairing we do another subset against dfMatchOverall, we can call this dfMatchOverallBest
dfMatchOverallBest <- subset(dfMatchOverall, dfMatchOverall$inGroupPairing %in% dfMatchOverallLineage1$inGroupPairing)

#To get the second lineage of each sister pairing we can subtract dfMatchOverallLineage1 by dfMatchOverall to get the bin uri's for the alternative pairings
dfMatchOverallLineage2 <- setdiff(dfMatchOverallBest$bin_uri, dfMatchOverallLineage1$bin_uri)
dfMatchOverallLineage2 <- as.data.frame(dfMatchOverallLineage2)
colnames(dfMatchOverallLineage2)[1] <- "bin_uri"
#merge with match overallbest to get all the relevant data for the second lineage
dfMatchOverallLineage2 <- merge(dfMatchOverallLineage2, dfMatchOverallBest, all.x = TRUE)
#Reorganization of columns for second lineage
dfMatchOverallLineage2 <- (dfMatchOverallLineage2[,c("inGroupPairing","record_id","bin_uri","inGroupDist","inGroupDistx1.3","medianLat","latMin","latMax","binSize","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","medianLon","indexNo")])

################
#Outgroup determination for each Pairing

#Now we can search for the best possible outgroupings for each pairing
#This will involve searching for outgroups relative to each pairing and finding one that is far enough away from both lineages

#First we can search for outgroupings relative to lineage 1
#We can use the indexNo to subset the dfGeneticDistanceStackDataframe according to indexes represented in lineage 1, this will be called dfBestOutGroupL1
#This will essentially limit our outgroup distances to those associated with lineage1
dfBestOutGroupL1 <- subset(dfGeneticDistanceStack, dfGeneticDistanceStack$ind %in% dfMatchOverallLineage1$indexNo)

#Then we can find which indices match lineage1 with bestoutgroup l1
outGroupCandidatesL1a <- foreach(i=1:nrow(dfMatchOverallLineage1)) %do% which(dfBestOutGroupL1$ind == dfMatchOverallLineage1$indexNo[i])
#Determining indices with correct outgroup distance 
outGroupCandidatesL1b <- foreach(i=1:nrow(dfMatchOverallLineage1)) %do% which(dfBestOutGroupL1$values >= dfMatchOverallLineage1$inGroupDistx1.3[i])
#Intersection of the two using mapply to find the correct outgroupings for each pairing
outGroupCandidatesL1c <- mapply(intersect,outGroupCandidatesL1a,outGroupCandidatesL1b)
#Unlist to make into one vector (less memory)
outGroupCandidatesL1c <- unlist(outGroupCandidatesL1c)
#Adding an rownum column to dfBestOutGroup, this represents the second index of each outgroup candidate
dfBestOutGroupL1$rownum<-seq.int(nrow(dfBestOutGroupL1))
#Then we can subset to rownum column based on the outgroup candidates
dfBestOutGroupL1 <- dfBestOutGroupL1[dfBestOutGroupL1$rownum %in% outGroupCandidatesL1c, ]

#Now for the second lineage, we go through the same process
dfBestOutGroupL2 <- subset(dfGeneticDistanceStack, dfGeneticDistanceStack$ind %in% dfMatchOverallLineage2$indexNo)
outGroupCandidatesL2a <- foreach(i=1:nrow(dfMatchOverallLineage2)) %do% which(dfBestOutGroupL2$ind == dfMatchOverallLineage2$indexNo[i])
outGroupCandidatesL2b <- foreach(i=1:nrow(dfMatchOverallLineage2)) %do% which(dfBestOutGroupL2$values >= dfMatchOverallLineage2$inGroupDistx1.3[i])
outGroupCandidatesL2c <- mapply(intersect,outGroupCandidatesL2a,outGroupCandidatesL2b)
outGroupCandidatesL2c <- unlist(outGroupCandidatesL2c)
dfBestOutGroupL2$rownum<-seq.int(nrow(dfBestOutGroupL2))
dfBestOutGroupL2 <- dfBestOutGroupL2[dfBestOutGroupL2$rownum %in% outGroupCandidatesL2c, ]

#Now we find the intersection between outGroupListL1 and outGroupListL2, this will determine which outgroup candidates are shared between both lineages
#Doing this for both L1 and L2
dfBestOutGroupL1 <- dfBestOutGroupL1[dfBestOutGroupL1$rownum %in% dfBestOutGroupL2$rownum, ]
#Then we break this dataframe down by index into a list
bestOutGroupListL1 <- lapply(unique(dfBestOutGroupL1$ind), function(x) dfBestOutGroupL1[dfBestOutGroupL1$ind == x,])
#We then find the value closest to the 1.3 divergence value
minOutGroupL1 <- sapply( bestOutGroupListL1 , function(x) min( x$values ) )
#make into a dataframe
dfMinOutGroupL1 <- as.data.frame(minOutGroupL1)
#Again subset to dfBestOutGroup to get indices of the closest outgrouping values
dfBestOutGroupL1 <- dfBestOutGroupL1[dfBestOutGroupL1$values %in% dfMinOutGroupL1$minOutGroupL1,]

#Again for L2
#For L2 its a bit different since we are choosing outgroups which would be shared between both lineages
#we simply subset dfBestOutGroupL2 rownum by dfBestOutGroupL1
dfBestOutGroupL2 <- dfBestOutGroupL2[dfBestOutGroupL2$rownum %in% dfBestOutGroupL1$rownum, ]

#Its possible to have more than one ideal outgrouping per pairing, so restricting to one outgrouping per pairing by setting head, n=1
#Once again doing this for both L1 and L2
dfBestOutGroupL1 <- by(dfBestOutGroupL1, dfBestOutGroupL1["ind"], head, n=1)
dfBestOutGroupL1 <- Reduce(rbind, dfBestOutGroupL1)
#Making sure that dfBestOutGroupOverall is ordered by rownum, this will be important for the next step
dfBestOutGroupL1 <- dfBestOutGroupL1[order(dfBestOutGroupL1$rownum),] 
#L2
dfBestOutGroupL2 <- by(dfBestOutGroupL2, dfBestOutGroupL2["ind"], head, n=1)
dfBestOutGroupL2 <- Reduce(rbind, dfBestOutGroupL2)
dfBestOutGroupL2 <- dfBestOutGroupL2[order(dfBestOutGroupL2$rownum),] 

#this will give the correct indexNo for the rownum column, doing this for both L1 and L2
j=0:(nrow(dfBestOutGroupL1)-1)
dfBestOutGroupL1$indexNo <- (dfBestOutGroupL1$rownum / (nrow(dfAllSeq))-j) * (nrow(dfAllSeq))
#then ordering by this indexNo, this is important for merging to the MatchOverall best dataframe
dfBestOutGroupL1<- dfBestOutGroupL1[order(dfBestOutGroupL1$indexNo),]
dfBestOutGroupL1 <- merge(dfBestOutGroupL1, dfAllSeq, by.x = "indexNo", by.y = "ind", all.x = TRUE)
#L2
j=0:(nrow(dfBestOutGroupL2)-1)
dfBestOutGroupL2$indexNo <- (dfBestOutGroupL2$rownum / (nrow(dfAllSeq))-j) * (nrow(dfAllSeq))
dfBestOutGroupL2<- dfBestOutGroupL2[order(dfBestOutGroupL2$indexNo),]
dfBestOutGroupL2 <- merge(dfBestOutGroupL2, dfAllSeq, by.x = "indexNo", by.y = "ind", all.x = TRUE) 

#Can rename certain columns to more closely resemble dfMatchOverallLineage dataframes except the outgroupdist column which would be unique of course
dfBestOutGroupL1 <- dfBestOutGroupL1[,c("bin_uri","record_id","values","medianLat","latMin","latMax","binSize","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","medianLon","indexNo")]
colnames(dfBestOutGroupL1)[3] <- "outGroupDist"
dfBestOutGroupL2 <- dfBestOutGroupL2[,c("bin_uri","record_id","values","medianLat","latMin","latMax","binSize","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","medianLon","indexNo")]
colnames(dfBestOutGroupL2)[3] <- "outGroupDist"

#Ordering values by order of values in minoutgroup dataframe to ensure we have the correct order
dfBestOutGroupL1 <- dfBestOutGroupL1[order(match(dfBestOutGroupL1[,3],dfMinOutGroupL1[,1])),]
#Then ordering outgroups for lineage 2 by ordering of lineage 1
dfBestOutGroupL2 <- dfBestOutGroupL2[order(match(dfBestOutGroupL2[,1],dfBestOutGroupL1[,1])),]

#adding an ingroup pairing column to each bestoutgroup dataframe
dfBestOutGroupL1$inGroupPairing <- dfMatchOverallLineage1$inGroupPairing
dfBestOutGroupL2$inGroupPairing <- dfMatchOverallLineage2$inGroupPairing

#Can now merge are outgroups to are associated lineages
#Each outgroup will be duplicated for each member of the pairing however the distances should be unqiue to represent each lineage

#.x beside headings are the ingroup data (except unique columns), scroll far enough to the right and headings with .y are the outgroup related columns
dfMatchOverallLineage1 <- merge(dfMatchOverallLineage1, dfBestOutGroupL1, by.x = "inGroupPairing", by.y = "inGroupPairing")
dfMatchOverallLineage2 <- merge(dfMatchOverallLineage2, dfBestOutGroupL2, by.x = "inGroupPairing", by.y = "inGroupPairing")

#Then we can revise dfMatchOverallBest to correctly reflect both pairings with correct outgroup distances
#We do this by using rbind to combine both lineages together into the dfMatchOverallBest dataframe
dfMatchOverallBest <-  rbind(dfMatchOverallLineage1, dfMatchOverallLineage2)
#Then order by inGroupPairing once again for a better organization of the dataframe
dfMatchOverallBest <- dfMatchOverallBest[order(dfMatchOverallBest$inGroupPairing),]

#As a last step for outgroup determination, if a suitable outgroup could not be found and is NA for a pairing, then we can filter these pairings out
noOutGroup<-which(is.na(dfMatchOverallBest$bin_uri.y))
if(length(noOutGroup) >0){
  dfMatchOverallBest<-dfMatchOverallBest[-noOutGroup,]}

################
#Ouput to CSV or TSV

#Should output the directory you have set in Rstudio

#Defining another variable to give a unique name to the CSV
#Name would be what you specify it to be, R will prompt you to insert a name for the file in the console:
#filename <- readline(prompt="")

#Once you have a name, uncomment one the of the write table commands and run that command to get a file output

#CSV
#write.table(dfMatchOverallBest, file=paste(filename, ".csv", sep=""), quote=FALSE, sep=',', col.names = NA)

#TSV
#write.table(dfMatchOverallBest, file=paste(filename, ".tsv", sep=""), quote=FALSE, sep='\t', col.names = NA)

################
