#################
#Authored by BIF IGP Group Tempura
#Major contributions by Matthew Orton to this script

#This program will allow for the generation of latitudinally separated sister pairings and associated outgroupings from ANY taxa and geographical region found on BOLD in one streamlined R pipeline!
#In this new and improved iteration, data is translated directly from a BOLD tsv file of the users choosing to a dataframe in R
#The generated sister pairs and outgroups can then be written to a csv or tsv and the file will appear in the current working directory of R

##################
#A few important tips:

#There are two options for parsing the tsv, you can either download the tsv directly from the BOLD API OR you can use a tsv you have previously downloaded

#Larger taxa, ex: aves take a several mins to process the tsv initially and run the sequence alignment. This can potentially consume a lot of working memory. I wouldnt run any taxa that is very large, (insecta for example) until we know what kind of memory resources it will consume.
#Ive found taxon=reptilia, geo=all to be really great for testing, it usually generates around 60 pairings and processing times are short

#Some tips for using the BOLD API since this is what is used to grab the relevant data we need: 
#To see details on how to use the bold API, go to http://www.boldsystems.org/index.php/resources/api?type=webservices
#Can add additional restrictions to url, for example &instituiton=Biodiversity Institute of Ontario|York University or &marker=COI-5P if you want to specifiy an institution or specific genetic marker in the url
#geo=all in the url means global but geo can be geo=Canada for example
#Can use | modifier in the url, for example &geo=Canada|Alaska would give data for both Canada and alaska, or taxon=Aves|Reptilia would yield give data for both Aves and Reptilia

#We recommend using RStudio: https://www.rstudio.com/home/, provides a nicer interface than just using R
#Large matrices and dataframes may not show all columns, you can overcome this by typing the command: utils::View() where you would insert the dataframe or matrix you want in the brackets of that command
#Pay attention to the seed you are using for the generation of randomly sampled bin members, keeping the seed the same will produce the same result each time, if you want random bin members generated each time however you can comment the seed out

#Important dataframes:
#dfMatchOverallBest is the finalized dataframe that contains all of the finalized pairings and outgroupings
#dfMatchOverall is the dataframe of ingroup pairings before bins with multiple pairings are eliminated based on distance
#dfMatchOverallPair is sort of an intermediate df between MatchOverall and MatchOverallBest that shows which pairing was chosen for each bin
#dfInitial is the dataframe first produced by the import from BOLD and is trimmed by lat, bin_uri etc.
#dfLatLon just contains relevant information for each bin: bin size, maximum lat, median lat, minimum lat, median lon
#dfBestOutGroup contains the associated outgroupings only but does have a column for which pairing its associated with
#dfLatitudeDistance and dfGeneticDistance are matrices converted into dataframes that show all possible distances between bins
#dfRandomBinSeq is the dataframe that contains information for all bin members that were randomly sampled from the sample function
#dfGeneticDistanceStack is dfGeneticDistance with all columns concatenated into one long column, it is used to grab index numbers for each pairing

#################
#Packages required

#For this we also need the foreach package for a few functions
install.packages("foreach")
library(foreach)
#For genetic distance, we use the ape package
install.packages("ape")
library(ape)
#Decided to use the data.tables package for one component of the code since I was having difficulties merging two df's
install.packages("data.table")
library(data.table)
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

#First we download the TSV and convert it into a dataframe, this URL is what is modified by the user and will determine the taxa, geographic region etc.
dfInitial <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Porifera&geo=all&format=tsv")

#If you want to run pre downloaded BOLD tsv's, this will let you choose a path to that tsv and parse
#tsvParseDoc <- file.choose()
#dfInitial <- read_tsv(tsvParseDoc)

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
containBin <- grep( "[:]", dfInitial$bin_uri)
dfInitial<-dfInitial[containBin,]

#Getting rid of any records that dont have sequence data (sometimes there are a few)
containNucleotides <- grep( "[ACGTN]", dfInitial$nucleotides)
dfInitial<-dfInitial[containNucleotides,]

#Modifying Bin column slightly to remove "BIN:"
dfInitial$bin_uri <- substr(dfInitial$bin_uri, 6 , 13)

#Dataframe Organization
dfInitial <- (dfInitial[,c("record_id","bin_uri","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","latNum","lonNum")])

#We first have to determine a median latitude for each bin so this will involve subdividing the dataframe into smaller dataframes by BIN

#First we can make a smaller dataframe with the columns we want for each bin - bin_uri, latnum, lonnum, record_id, if we didnt do this, the binList would consume a huge amount of memory
dfBinList <- (dfInitial[,c("record_id","bin_uri","latNum","lonNum")])
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

#We also will randomly select sequences from each bin of dfInitial with all of its sequence and taxonomic data, this is because for larger taxa, it may become too computionally intensive to run all sequences of each bin
#sapply will iterate through each bin and sample one of its members from the initial dataframe randomly
#Had set a seed for testing - meaning it will produce the same result each time, if you want a random result each time, comment seed out
set.seed(15);
record_id <- sapply( binList , function(x) sample( x$record_id, 1, replace = FALSE ) )
dfRandomBinSeq <- data.frame(record_id)
dfRandomBinSeq <- merge(dfRandomBinSeq, dfInitial, by.x = "record_id", by.y = "record_id")

#Can now merge dfLatLon with RandomBinseq (but retain the name RandomBinSeq) to create one dataframe containg all pertinent data we need and can reference back to from the match overall dataframe produced later on
dfRandomBinSeq <- merge(dfRandomBinSeq, dfLatLon, by.x = "bin_uri")
#Getting rid of latNum and lonNum since we dont need these anymore, we only need the median Lat/Lon values 
dfRandomBinSeq <- (dfRandomBinSeq[,c("bin_uri","binSize","record_id","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","medianLat","latMin","latMax","medianLon")])
#Adding an index column to reference later with Match overall dataframe
dfRandomBinSeq$ind <- row.names(dfRandomBinSeq)

#Can then determine latitude differences, originally put into a matrix but transformed into a df
#Can be found easily with dist function and put into a matrix
distLat <- dist(dfRandomBinSeq$medianLat)
matrixLatitudeDistance <- as.matrix( dist(dfRandomBinSeq$medianLat) )
#Convert to dataframe since its easier to manipulate, this dataframe will be used further down
dfLatitudeDistance <-as.data.frame(matrixLatitudeDistance)

#Converting sequences to DNAStringSet format
dnaStringSet <- DNAStringSet(dfRandomBinSeq$nucleotides)

#Multiple sequence alignment using msa package - 3 different algorithms can be used for this: ClustalW, ClustalOmega, MUSCLE, using the default of ClustalW for testing
#Sequence alignment package can be found here: http://master.bioconductor.org/packages/3.2/bioc/vignettes/msa/inst/doc/msa.pdf
#This could take several minutes to even hours in some cases depending on the taxa
seqAlignment <- msa(dnaStringSet)

#Conversion to DNAbin format before genetic distance computation matrix
DNAbin <- as.DNAbin(seqAlignment)

#Now for the computation of genetic distance, several models can be used - "raw", "N", "TS", "TV", "JC69", "K80" (the default), "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet"
#Details on each model can be found here: http://www.inside-r.org/packages/cran/ape/docs/dist.dna
#Starting with the K80 model 
#Note that this has only been tested with K80 currently so unsure as to how it will perform with other models
matrixGeneticDistance <- dist.dna(DNAbin, model = "K80", as.matrix = TRUE, pairwise.deletion = TRUE)
#convert to dataframe
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
dfMatchOverall <- (dfMatchOverall[,c("inGroupPairing","record_id","bin_uri","values","inGroupDistx1.3","medianLat","latMin","latMax","binSize","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","medianLon","ind")])
colnames(dfMatchOverall)[4] <- "inGroupDist"
colnames(dfMatchOverall)[26] <- "indexNo"

#For bins that have multiple pairings we can select which pairing has the smallest divergence to its other ingroup member since the pairings are ordered according to distance
#Basically we are taking the best possible pairing per bin based on divergence
#This dataframe will be dfMatchOverallPair

#This command selects the top pairing for each bin since the bins are ordered already
dfMatchOverallPair <- by(dfMatchOverall, dfMatchOverall["bin_uri"], head, n=1)
#Makes another dataframe
dfMatchOverallPair <- Reduce(rbind, dfMatchOverallPair)
dfMatchOverallPair <- subset(dfMatchOverallPair, duplicated(dfMatchOverallPair$inGroupPairing))
#Can then take our best pairings and subset them against the original dfMatchOverall dataframe to get the finalized dataframe of ingroup pairings which will be called dfMatchOverallBest
dfMatchOverallBest <- subset(dfMatchOverall, dfMatchOverall$inGroupPairing %in% dfMatchOverallPair$inGroupPairing)

#We now have the best possible pairings from each bin based on both latitude and distance criteria!

#We now have to determine the best possible outgroupings

#Can find suitable outgroupings based on the calculated mimumum genetic distances of the ingroupdist1.3x column and indexNo columns
#First we can use the indexNo to subset the dfGeneticDistanceStackDataframe according to indexes represented in the final pairings, this will be called dfBestOutGroup
#This will essentially limit our outgroup bin distances to those associated with the ingroup bin distances
dfBestOutGroup <- subset(dfGeneticDistanceStack, dfGeneticDistanceStack$ind %in% dfMatchOverallPair$indexNo)

#Determining shared indices between pairings and dfBestOutGroup
outGroupCandidates <- foreach(i=1:nrow(dfMatchOverallPair)) %do% which(dfBestOutGroup$ind == dfMatchOverallPair$indexNo[i])
#Determining indices with correct outgroup distance 
outGroupCandidates1 <- foreach(i=1:nrow(dfMatchOverallPair)) %do% which(dfBestOutGroup$values >= dfMatchOverallPair$inGroupDistx1.3[i])
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
colnames(dfRandomBinSeq)[23] <- "identity"
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
dfBestOutGroup <- dfBestOutGroup[,c("bin_uri","record_id","values","medianLat","latMin","latMax","binSize","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","medianLon","identity")]
colnames(dfBestOutGroup)[3] <- "outGroupDist"
colnames(dfBestOutGroup)[24] <- "indexNo"

#Now just need to associate our outgroup with the pairings by adding another column to the outgroup df 
#Just ordering values by order of values in ClosestOutGroup to ensure we have the correct order
dfBestOutGroup <- dfBestOutGroup[order(match(dfBestOutGroup[,3],dfClosestOutGroup[,1])),]

#Appending the ingrouppairing column to the Outgroup dataframe so we know which pairing its associated with
dfBestOutGroup$associatedInGroup <- dfMatchOverallPair$inGroupPairing

#Can now merge eveything together into the MatchOverallBest, this will act as our finalized dataframe
#Note that each outgroup will be duplicated for each member of the pairing

#.x beside headings are the ingroup data (except unique columns), scroll far enough to the right and headings with .y are the outgroup related columns
dfMatchOverallBest <- merge(dfMatchOverallBest, dfBestOutGroup, by.x = "inGroupPairing", by.y = "associatedInGroup")

#Note there sometimes will be a few pairings that dont have a suitable outgroup(couldnt find a bin with the suitable outgroup distance) so we can filter those out
noOutGroupFilter<-which(is.na(dfMatchOverallBest$bin_uri.y))
if(length(noOutGroupFilter) >0){
  dfMatchOverallBest<-dfMatchOverallBest[-noOutGroupFilter,]}

#As one last step, we have to determine the outgroup distance to the other ingroup lineage, since we only have a distance for one lineage
#we determine the bin_uri's of the other lineage for each pairing
oppositeLineage <- setdiff(dfMatchOverallBest$bin_uri.x, dfMatchOverallPair$bin_uri)
#Make into a dataframe
dfOppositeLineage <- as.data.frame(oppositeLineage)
#Merge this with MatchoverallBest to get all the relevant data we need
dfOppositeLineage <- merge(dfMatchOverallBest, dfOppositeLineage, by.x = "bin_uri.x", by.y = "oppositeLineage")
#Subset according to our Gentic distance matrix to grab the distance for each lineage
outGroupDistOppositeLineage <- foreach(k=1:nrow(dfOppositeLineage)) %do% dfGeneticDistance[dfOppositeLineage$indexNo.x[k],dfOppositeLineage$indexNo.y[k]]
#Add it to the dfOppositeLineage df
dfOppositeLineage$outGroupDist <- outGroupDistOppositeLineage
#Add outgroup to first lineage
dfMatchOverallPair <- merge(dfMatchOverallPair, dfBestOutGroup, by.x = "inGroupPairing", by.y = "associatedInGroup")
#Add second lineage to first lineage and revise MatchOverallBest
dfMatchOverallBest <-  rbind(dfMatchOverallPair, dfOppositeLineage)
#Order by pairing once again
dfMatchOverallBest<- dfMatchOverallBest[order(dfMatchOverallBest$inGroupPairing),] 

#Now we have both outgroups and suitable pairings with distances from outgroup for both lineages of a pairing!

############
#Ouput to CSV or TSV

#Should output the directory you have set in Rstudio

#Defining another variable to give a unqiue name to the CSV
#Name would be what you specify it to be, R will prompt you to insert a name for the file in the console:
#filename <- readline(prompt="")

#Once you have a name, uncomment one the of the write table commands and run that command to get a file output

#CSV
#write.table(dfMatchOverallBest, file=paste(filename, ".csv", sep=""), quote=FALSE, sep=',', col.names = NA)

#TSV
#write.table(dfMatchOverallBest, file=paste(filename, ".tsv", sep=""), quote=FALSE, sep='\t', col.names = NA)
#############
