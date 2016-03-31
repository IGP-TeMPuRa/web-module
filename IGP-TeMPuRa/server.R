# server.R

#Packages required
library(shiny)
library(leaflet)
library(data.table)
library(DT)
#We need the foreach package for several functions that require iteration over dataframe rows
#install.packages("foreach")
library(foreach)
#For genetic distance determination using the TN93 model, we use the ape package
#install.packages("ape")
library(ape)
#Speeds up parsing of the tsv file with read_tsv function
#install.packages("readr")
library(readr)
#For sequence alignments we need the biostrings (DNAStringSet function) and msa packages, run each of these commands individually, sometimes it skips the libraries
source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
library("Biostrings")
#biocLite("msa")
library("msa")
#For overlapping latitude regions we need the Desctools package
#install.packages("DescTools")
library(DescTools)
#Also adding data tables for table merging
#install.packages("data.table")
library(data.table)
#For plotting of relative outgroup distances between lineages we will also need ggplot2
require(ggplot2)

shinyServer(function(input, output, session) {
  
  #Defining variables outside of the reactive so that it can be stored temporarily when the application is running
  var1 <- NULL
  var2 <- NULL
  var3 <- NULL
  URL <- NULL
  dfInitial <- NULL
  containLat <- NULL
  latNum <- NULL
  lonNum <- NULL
  containBin <- NULL
  containNucleotides <- NULL
  dfBinList <- NULL
  binList <- NULL
  medianLat <- NULL
  medianLon <- NULL
  latMin <- NULL
  latMax <- NULL
  binSize <- NULL
  dfLatLon <- NULL
  largeBin <- NULL
  dfConsensus <- NULL
  binNumberConsensus <- NULL
  dfNonConsensus <- NULL
  largeBinList <- NULL
  dnaStringSet1 <- NULL
  alignment1 <- NULL
  consensusSeq <- NULL
  dfConsensus <- NULL
  dfAllSeq <- NULL
  distLat <- NULL
  matrixLatitudeDistance <- NULL
  dfLatitudeDistance <- NULL
  alignmentSequences <- NULL
  dnaStringSet2 <- NULL
  alignment2 <- NULL
  DNAbin <- NULL
  matrixGeneticDistance <- NULL
  dfGeneticDistance <- NULL
  dfGeneticDistanceStack <- NULL
  geneticDistanceMatchI <- NULL
  latitudeDistanceMatch <- NULL
  matchOverall <- NULL
  dfMatchOverall <- NULL
  dfMatchOverallLineage1 <- NULL
  dfMatchOverallBest <- NULL
  dfMatchOverallLineage2 <- NULL
  rangeThreshold <- NULL
  rangeL1 <- NULL
  rangeL2 <- NULL
  overlapValue <- NULL
  rangeOverlapCheck <- NULL
  overlapInd <- NULL
  dfDistancePair <- NULL
  closerBinRow <- NULL
  closerBinColumn <- NULL
  dfBestOutGroupL1 <- NULL
  outGroupCandidatesL1a <- NULL
  outGroupCandidatesL1b <- NULL
  outGroupCandidatesL1c <- NULL
  dfBestOutGroupL2 <- NULL
  outGroupCandidatesL2a <- NULL
  outGroupCandidatesL2b <- NULL
  outGroupCandidatesL2c <- NULL
  indNum <- NULL
  noOutGroup <- NULL
  
  urlInput <- reactive({
    var1 <<- "http://www.boldsystems.org/index.php/API_Public/combined?taxon="
    var2 <<- "&geo="
    var3 <<- "&format=tsv"
    URL <<- paste0(var1, input$taxonomy, var2, input$geography, var3)
    dfInitial <<- read_tsv(URL)
    #Dataframe Filtering and Reorganization
    
    #Filtering this df according to the relevant columns we need
    dfInitial <<- (dfInitial[,c("recordID","bin_uri","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","lat","lon","nucleotides")])
    colnames(dfInitial)[1] <<- "record_id"
    
    #Removing sequences with no latitude values, filtering according to lat since we only really need lat for the analysis
    containLat <<- grep( "[0-9]", dfInitial$lat)
    dfInitial<<-dfInitial[containLat,]
    
    #Next we have to convert lat column to num instead of chr type, this will become important later on for median latitude determination
    latNum <<- with(dfInitial, as.numeric(as.character(lat))) 
    dfInitial$latNum <<- latNum
    
    #Can do lon as well to get numeric values instead of characters
    lonNum <<- with(dfInitial, as.numeric(as.character(lon))) 
    dfInitial$lonNum <<- lonNum
    
    #First identifying missing bins and eliminating rows with missing bin_uri's since bin is a big indicator of sequence quality
    #Grep by colon since every record with a bin identifier will have this
    containBin <<- grep( "[:]", dfInitial$bin_uri)
    dfInitial<<-dfInitial[containBin,]
    
    #Getting rid of any records that dont have sequence data (sometimes there are a few)
    containNucleotides <<- grep( "[ACGTN]", dfInitial$nucleotides)
    dfInitial<<-dfInitial[containNucleotides,]
    
    #Modifying Bin column slightly to remove "BIN:"
    dfInitial$bin_uri <<- substr(dfInitial$bin_uri, 6 , 13)
    
    #Dataframe Reorganization
    dfInitial <<- (dfInitial[,c("record_id","bin_uri","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","latNum","lonNum")])
    
    #Bin Stats and Median Latitude/Longitude Determination per bin
    
    #First we can make a smaller dataframe with the columns we want for each bin - bin_uri, latnum, lonnum, record_id, if we didnt do this, the binList would consume a huge amount of memory
    dfBinList <<- (dfInitial[,c("record_id","bin_uri","latNum","lonNum","nucleotides")])
    #Create groupings by bin with each grouping representing a different bin_uri
    #Each element of this list represents a bin with subelements representing the various columns of the initial dataframe created and the information is grouped by bin 
    binList <<- lapply(unique(dfBinList$bin_uri), function(x) dfBinList[dfBinList$bin_uri == x,])
    
    #Now to determine a median latitude for each bin
    medianLat <<- sapply( binList , function(x) median( x$latNum ) )
    
    #We also need a median longitude for each if we are going to plot on a map for a visual interface
    medianLon <<- sapply( binList , function(x) median( x$lonNum ) )
    
    #we can also take a few other important pieces of data regarding each bin using sapply including number of record_ids to a bin and latitudinal min and max of each bin
    latMin <<- sapply( binList , function(x) min( x$latNum ) )
    latMax <<- sapply( binList , function(x) max( x$latNum ) )
    binSize <<- sapply( binList , function (x) length( x$record_id ) )
    
    #Dataframe of our median lat values, this will be used in our final dataframe
    dfLatLon <<- data.frame(medianLat)
    
    #Adding bin_uri, median longitude, latMin, latMax and binSize to dataframe with medianLat
    dfLatLon$bin_uri <<- c(unique(dfInitial$bin_uri))
    dfLatLon$medianLon <<- c(medianLon)
    dfLatLon$latMin <<- c(latMin)
    dfLatLon$latMax <<- c(latMax)
    dfLatLon$binSize <<- c(binSize)
    
    #Can also convert to 180 degree latitude scale, convert longitude as well
    dfLatLon$medianLat <<- dfLatLon$medianLat + 90
    dfLatLon$latMax <<- dfLatLon$latMax + 90
    dfLatLon$latMin <<- dfLatLon$latMin + 90
    dfLatLon$medianLon <<- dfLatLon$medianLon + 90
    
    #Merging LatLon to BinList for the sequence alignment step
    dfBinList <<- merge(dfBinList, dfLatLon, by.x = "bin_uri", by.y = "bin_uri")
    
    ###############
    #Sequence alignments and consensus sequences for each bin with bins larger than 1 member
    
    #First we have to subset dfBinList to find bins with more than one member since these bins will need consensus sequences
    largeBin <<-which(dfBinList$binSize > 1)
    #If there is at least one bin with more than one member, then a dataframe dfConsensus will be created with those bins
    if(length(largeBin) >0){
      dfConsensus <<- dfBinList[largeBin,]
      
      #Also need to find the number of unique bins in dfConsensus
      binNumberConsensus <<- unique(dfConsensus$bin_uri)
      binNumberConsensus <<- length(binNumberConsensus)
      
      #We also have to create another separate dataframe with bins that only have one member
      dfNonConsensus <<- dfBinList[-largeBin,]
      
      #We then take the dfConsensus sequences and break it down into a list with each element being a unique bin
      largeBinList <<- lapply(unique(dfConsensus$bin_uri), function(x) dfConsensus[dfConsensus$bin_uri == x,])
      
      #Convert all of the sequences in this list to dnaStringSet format for the alignment step
      dnaStringSet1 <<- sapply( largeBinList, function(x) DNAStringSet(x$nucleotides) )
      
      #Multiple sequence alignment using msa package - 3 different algorithms can be used for this: ClustalW, ClustalOmega, MUSCLE, using the default of ClustalW for testing
      #Sequence alignment package can be found here: http://master.bioconductor.org/packages/3.2/bioc/vignettes/msa/inst/doc/msa.pdf
      #Using the default parameters of the msa package (ClustalW)
      #Run a multiple sequence alignment on each element of the dnaStringSet1 list, can achieve this with the foreach package
      alignment1 <<- foreach(i=1:binNumberConsensus) %do% msa(dnaStringSet1[[i]])
      
      #We can then extract the consensus sequence from each individual alignment
      consensusSeq <<- foreach(i=1:binNumberConsensus) %do% consensusString(alignment1[[i]])
      
      #Adding our consensusSeq to dfConsensus 
      dfConsensus <<- by(dfConsensus, dfConsensus["bin_uri"], head, n=1)
      dfConsensus <<- Reduce(rbind, dfConsensus)
      dfConsensus$nucleotides <- consensusSeq
      #Labeling all consensus sequences as such in the record_id since they would not have a specific record_id
      dfConsensus$record_id <<- "consensus sequence"
      
      #Append our nonconsensus to our consensus, will make a new dataframe for this - dfAllSeq
      #now we have all of the sequences we need for the next alignment of all sequences
      dfAllSeq <<-  rbind(dfConsensus, dfNonConsensus)
      #Can then merge this with dfInitial to get all of the relevant data we need
      #merging to dfInitial to gain all the relevant taxanomic data
      #of course consensus sequences will not have any specific taxonomic data but will still retain their bin_uri for identification
      dfAllSeq <<- merge(dfAllSeq, dfInitial, all.x = TRUE)
      #Renaming and reorganizing the dataframe
      dfAllSeq <<- (dfAllSeq[,c("bin_uri","binSize","record_id","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","medianLat","latMin","latMax","medianLon")])
      #Adding an index column to reference later with Match overall dataframe
      dfAllSeq$ind <<- row.names(dfAllSeq)
      
    } else {
      #Else if there are no bins with more than one member than we would simply merge latlon with intial to get dfAllSeq and thus all of the sequences we want
      dfAllSeq <<- merge(dfLatLon, dfInitial, by.x = "record_id", by.y = "record_id")
      dfAllSeq <<- (dfAllSeq[,c("bin_uri","binSize","record_id","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","medianLat","latMin","latMax","medianLon")])
      dfAllSeq$ind <<- row.names(dfAllSeq)
    }
    
    #############
    #Latitude Distance Determination
    
    #Can be determined easily with the dist function which will determine latitudinal differences between all bins
    distLat <<- dist(dfAllSeq$medianLat)
    #Then we convert this to a matrix
    matrixLatitudeDistance <<- as.matrix( dist(dfAllSeq$medianLat) )
    #Then convert to dataframe since its easier to manipulate, this dataframe will be used further down (in testing I found matrices tend to be less reliable for further manipulations)
    dfLatitudeDistance <<-as.data.frame(matrixLatitudeDistance)
  })

  
  
  distanceModels <- reactive({
    urlInput()
    #############
    #Multiple Sequence Alignment of All Sequences and Pairwise Distance determination with TN93
    
    #Lets first start off by identifying our reference sequence
    #we can make a smaller dataframe with the name of the taxa as one column and the sequence as another column
    #Can call this dfRefSeq
    
    #We also have to convert to type character for this to work
    alignmentSequences <<- as.character(dfAllSeq$nucleotides)
    
    #Converting all sequences in dfAllSeq to DNAStringSet format, this is the format required for the alignment
    dnaStringSet2 <<- DNAStringSet(alignmentSequences)
    
    #Run a multiple sequence alignment of all sequences both consensus and nonconsensus
    #As mentioned using the default of ClustalW
    #This could take several minutes to even hours in some cases depending on the taxa
    alignment2 <<- msa(dnaStringSet2)
    
    #Conversion to DNAbin format before genetic distance matrix
    DNAbin <<- as.DNAbin(alignment2)
    
    #Now for the computation of genetic distance, several models can be used - "raw", "N", "TS", "TV", "JC69", "K80" (the default), "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet"
    #Details on each model can be found here: http://www.inside-r.org/packages/cran/ape/docs/dist.dna
    #Using the TN93 model for our data
    matrixGeneticDistance <<- dist.dna(DNAbin, model = input$distanceModels, as.matrix = TRUE, pairwise.deletion = TRUE)
    #convert to dataframe
    dfGeneticDistance <<-as.data.frame(matrixGeneticDistance)
    #Putting it into a stack (each column concatenated into one long column of indexes and values) so it can be easily subsetted
    dfGeneticDistanceStack <<-stack(dfGeneticDistance)
  })
  
  filterInput <- reactive({
    distanceModels()
    ##############
    #Finding appropriate pairings according to latitude and distance criteria

    #So now we have two dataframes, we can find ideal matchings based on less than 15% divergence and 20 degrees latitude separation
    #These values can easily be edited to add more or less stringency to the matches
    #Will produce lists with indexes of each match according to our set criteria
    geneticDistanceMatchI <<- which(dfGeneticDistance<= input$genetic)
    latitudeDistanceMatch <<- which(dfLatitudeDistance>= input$latitude)

    #Match the dfs against each other to find the sister pairs using intersect
    matchOverall <<- intersect(geneticDistanceMatchI,latitudeDistanceMatch)

    #Use these index matches and reference against the genetic distance stack dataframe to subset it according to matches
    dfMatchOverall <<- dfGeneticDistanceStack[c(matchOverall), ]
    #Each duplication of genetic distance values corresponds to a pairing
    #Will now merge the randombinseq dataframe to the pairings generated
    dfMatchOverall <<- merge(dfAllSeq, dfMatchOverall, by.x = "ind", by.y = "ind")

    #Pairings will be ordered according to genetic distance, least divergent pairing to most divergent pairing
    dfMatchOverall <<- dfMatchOverall[order(dfMatchOverall$values),]
  })
  
  outgroupInput <- reactive({
    filterInput()
    #Then we can multiply these distance values by 1.3 to determine the minimum outgroup distance from the pairings, we put this in another column in the matchOverall dataframe
    #This minimum outgroup distance could also be a user adjustable parameter
    dfMatchOverall$inGroupDistx1.3 <<- dfMatchOverall$values * input$outgroups

    #Also grouping dataframe every 2 rows to reflect each unique pairing/matching, pairing column will give a number value for each pairing ordered
    dfMatchOverall$inGroupPairing <<- rep(1:(nrow(dfMatchOverall)/2), each = 2)

    #Reorganizing and renaming some columns in MatchOverall dataframe to make more easily readable
    dfMatchOverall <<- (dfMatchOverall[,c("inGroupPairing","record_id","bin_uri","values","inGroupDistx1.3","medianLat","latMin","latMax","binSize","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","medianLon","ind")])
    colnames(dfMatchOverall)[4] <<- "inGroupDist"
    colnames(dfMatchOverall)[26] <<- "indexNo"

    ##############
    #Taking the Best Possible Pairings (no bin duplicated in any pairing) and Creating Dataframes for each Pairing Lineage and each Complete Pairing (both lineages)

    #For bins that have multiple pairings associated with them we can select which pairing has the smallest divergence to its other ingroup member since the pairings are ordered according to distance
    #Basically we are taking the best possible pairing per bin based on divergence such that a bin would never be duplicated in multiple pairings
    #This dataframe will be dfMatchOverallLineage1
    dfMatchOverallLineage1 <<- by(dfMatchOverall, dfMatchOverall["bin_uri"], head, n=1)
    #This dataframe represents the first lineage of each sister pairing
    dfMatchOverallLineage1 <<- Reduce(rbind, dfMatchOverallLineage1)
    dfMatchOverallLineage1 <<- subset(dfMatchOverallLineage1, duplicated(dfMatchOverallLineage1$inGroupPairing))

    #To get both lineages for each pairing we do another subset against dfMatchOverall, we can call this dfMatchOverallBest
    dfMatchOverallBest <<- subset(dfMatchOverall, dfMatchOverall$inGroupPairing %in% dfMatchOverallLineage1$inGroupPairing)

    #To get the second lineage of each sister pairing we can subtract dfMatchOverallLineage1 by dfMatchOverall to get the bin uri's for the alternative lineage
    dfMatchOverallLineage2 <<- setdiff(dfMatchOverallBest$bin_uri, dfMatchOverallLineage1$bin_uri)
    dfMatchOverallLineage2 <<- as.data.frame(dfMatchOverallLineage2)
    colnames(dfMatchOverallLineage2)[1] <<- "bin_uri"
    #merge with match overallbest to get all the relevant data for the second lineage
    dfMatchOverallLineage2 <<- merge(dfMatchOverallLineage2, dfMatchOverallBest, all.x = TRUE)
    #Reorganization of columns for second lineage
    dfMatchOverallLineage2 <<- (dfMatchOverallLineage2[,c("inGroupPairing","record_id","bin_uri","inGroupDist","inGroupDistx1.3","medianLat","latMin","latMax","binSize","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","medianLon","indexNo")])

    ##############
    #Eliminating Pairings based on Overlapping Latitudinal Range

    #Next we can work on establishing latitudinal ranges for each pairing
    #If the two lineages of a pairing have overlapping latitude regions of greater than 25% then we would not consider that pairing as a viable pairing
    #This is because we want each lineage of a pairing to meet an appropriate difference in latitude

    #we can define the overlap range threshold as 25% of the latitude range of L1
    #If an overlap is greater than this value we would discard with this pairing
    #Of course this value could be easily modified to add more or less stringency to the script
    rangeThreshold <<- foreach(l=1:nrow(dfMatchOverallLineage1)) %do% ((dfMatchOverallLineage1$latMax[l] - dfMatchOverallLineage1$latMin[l]) * 0.25)

    #Define our latitude ranges for each lineage of a pairing
    rangeL1 <<- foreach(l=1:nrow(dfMatchOverallLineage1)) %do% range(dfMatchOverallLineage1$latMax[l], dfMatchOverallLineage1$latMin[l])
    rangeL2 <<- foreach(l=1:nrow(dfMatchOverallLineage1)) %do% range(dfMatchOverallLineage2$latMax[l], dfMatchOverallLineage2$latMin[l])

    #Then we can determine the overlap region between them using the Overlap function from the Desctools package
    #Overlap will return an absolute value so we dont have to worry about negatives for overlap values
    overlapValue <<- foreach(l=1:nrow(dfMatchOverallLineage1)) %do% Overlap(rangeL1[[l]], rangeL2[[l]])

    #Then if there is a range overlap between two lineages in a pairing, we can determine if this overlap is actually larger than the 25% value of rangeThreshold for each individual pairing
    rangeOverlapCheck <<- foreach(l=1:nrow(dfMatchOverallLineage1)) %do% which(rangeThreshold[[l]]<overlapValue[[l]])

    #Then overlaps values exceeding that 25% value we set should be returned as an integer of 1, if an overlap does not exceed this value, then it will return a value of 0
    #If there is a value that meets these criteria, its associated pairing will be removed from the MatchOverallBest and lineage dataframes (we will retain it in dfMatchOverall for reference)
    #This will name each element of rangeOverlapCheck with the inGroupPairing number to identify the pairing we need to eliminate
    names(rangeOverlapCheck) <<- paste0(dfMatchOverallLineage1$inGroupPairing)
    #Identify which pairing in rangeOverlapCheck is greater than 0
    overlapInd <<- which(rangeOverlapCheck>0)
    #if overlapInd is not empty:
    if(length(overlapInd)>0){
      #Eliminate based on that pairing number(s) for MatchOverallLineage1
      dfMatchOverallLineage1 <<- dfMatchOverallLineage1[-overlapInd,]
      #Now subset to both dfMatchOverallBest and dfMatchOverallLineage2
      dfMatchOverallBest <<- subset(dfMatchOverallBest, dfMatchOverallBest$inGroupPairing %in% dfMatchOverallLineage1$inGroupPairing)
      dfMatchOverallLineage2 <<- subset(dfMatchOverallLineage2, dfMatchOverallLineage2$inGroupPairing %in% dfMatchOverallLineage1$inGroupPairing)
    }


    ###############
    #Identifying and Averaging PseudoReplicates

    #To check for the phylogenetics problem of pseduoreplication we can generate another smaller distance matrix with our selected pairings only
    #If a bin from one pairing is actually closer to a bin from another pairing (as opposed to its paired lineage) than we would have to average the results of those two pairings

    #First lets order our lineage according to distance (makes things easier to understand and ensure we have the right ordering)
    dfMatchOverallLineage1 <<- dfMatchOverallLineage1[order(dfMatchOverallLineage1$inGroupDist),]
    dfMatchOverallLineage2 <<- dfMatchOverallLineage2[order(dfMatchOverallLineage2$inGroupDist),]

    #To do this we can subset our original pairwise distance matrix with the pairing indices of each lineage only
    #We will call this new dataframe dfDistancePair
    dfDistancePair <<- dfGeneticDistance[dfMatchOverallLineage1$indexNo,dfMatchOverallLineage2$indexNo]

    #For each column and row of this matrix, if a distance value is lower than the pairwise distances of each pairing than we know there is another bin which is actually closer
    #We check for this by using our pairwise distances from the dfMatchOverallLineages dataframes
    #First we check if another bin is closer for each row of the distance pair matrix
    closerBinRow <<- foreach(j=1:nrow(dfDistancePair)) %do% which(dfDistancePair[dfMatchOverallLineage1$indexNo,dfMatchOverallLineage2$indexNo[j]]<dfMatchOverallLineage1$inGroupDist[j])

    #Then we check if another bin is closer for each column of the distance pair matrix
    closerBinColumn <<- foreach(j=1:nrow(dfDistancePair)) %do% which(dfDistancePair[dfMatchOverallLineage1$indexNo[j],dfMatchOverallLineage2$indexNo]<dfMatchOverallLineage1$inGroupDist[j])

    #We then name each element of closerBinRow with the inGroupPairing number to identify the pairings we need to average
    names(closerBinRow) <<- paste0(dfMatchOverallLineage1$inGroupPairing)

    #we then name each element of closerBinColumn with the inGroupPairing number to identify the pairings we need to average
    names(closerBinColumn) <<- paste0(dfMatchOverallLineage1$inGroupPairing)

    ################
    #Outgroup determination for each Pairing

    #Now we can search for the best possible outgroupings for each pairing
    #This will involve searching for outgroups relative to each pairing and finding one that is far enough away from both lineages

    #First we can search for outgroupings relative to lineage 1
    #We can use the indexNo to subset the dfGeneticDistanceStackDataframe according to indexes represented in lineage 1, this will be called dfBestOutGroupL1
    #This will essentially limit our outgroup distances to those associated with lineage1
    dfBestOutGroupL1 <<- subset(dfGeneticDistanceStack, dfGeneticDistanceStack$ind %in% dfMatchOverallLineage1$indexNo)
    #We also order by Lineage1
    dfBestOutGroupL1 <<- dfBestOutGroupL1[order(match(dfBestOutGroupL1[,2],dfMatchOverallLineage1[,26])),]

    #Then we can find which indices match lineage1 with bestoutgroup l1
    outGroupCandidatesL1a <<- foreach(i=1:nrow(dfMatchOverallLineage1)) %do% which(dfBestOutGroupL1$ind == dfMatchOverallLineage1$indexNo[i])
    #Determining indices with correct outgroup distance
    outGroupCandidatesL1b <<- foreach(i=1:nrow(dfMatchOverallLineage1)) %do% which(dfBestOutGroupL1$values >= dfMatchOverallLineage1$inGroupDistx1.3[i])
    #Intersection of the two using mapply to find the correct outgroupings for each pairing
    outGroupCandidatesL1c <<- mapply(intersect,outGroupCandidatesL1a,outGroupCandidatesL1b)
    #Unlist to make into one vector (less memory)
    outGroupCandidatesL1c <<- unlist(outGroupCandidatesL1c)
    #Adding an rownum column to dfBestOutGroup, this represents the second index of each outgroup candidate
    dfBestOutGroupL1$rownum<<-seq.int(nrow(dfBestOutGroupL1))
    #Then we can subset to rownum column based on the outgroup candidates
    dfBestOutGroupL1 <<- dfBestOutGroupL1[dfBestOutGroupL1$rownum %in% outGroupCandidatesL1c, ]

    #Now for the second lineage, we go through the same process
    dfBestOutGroupL2 <<- subset(dfGeneticDistanceStack, dfGeneticDistanceStack$ind %in% dfMatchOverallLineage2$indexNo)
    dfBestOutGroupL2 <<- dfBestOutGroupL2[order(match(dfBestOutGroupL2[,2],dfMatchOverallLineage2[,26])),]
    outGroupCandidatesL2a <<- foreach(i=1:nrow(dfMatchOverallLineage2)) %do% which(dfBestOutGroupL2$ind == dfMatchOverallLineage2$indexNo[i])
    outGroupCandidatesL2b <<- foreach(i=1:nrow(dfMatchOverallLineage2)) %do% which(dfBestOutGroupL2$values >= dfMatchOverallLineage2$inGroupDistx1.3[i])
    outGroupCandidatesL2c <<- mapply(intersect,outGroupCandidatesL2a,outGroupCandidatesL2b)
    outGroupCandidatesL2c <<- unlist(outGroupCandidatesL2c)
    dfBestOutGroupL2$rownum<<-seq.int(nrow(dfBestOutGroupL2))
    dfBestOutGroupL2 <<- dfBestOutGroupL2[dfBestOutGroupL2$rownum %in% outGroupCandidatesL2c, ]

    #Now we find the intersection between outGroupListL1 and outGroupListL2, this will determine which outgroup candidates are shared between both lineages
    dfBestOutGroupL1 <<- dfBestOutGroupL1[dfBestOutGroupL1$rownum %in% dfBestOutGroupL2$rownum, ]
    dfBestOutGroupL2 <<- dfBestOutGroupL2[dfBestOutGroupL2$rownum %in% dfBestOutGroupL1$rownum, ]

    #Its possible to have more than one outgrouping per pairing that meets the divergence criteria, so restricting to one outgrouping per pairing by setting head, n=1
    #Once again doing this for both L1 and L2
    dfBestOutGroupL1 <<- by(dfBestOutGroupL1, dfBestOutGroupL1["ind"], head, n=1)
    dfBestOutGroupL1 <<- Reduce(rbind, dfBestOutGroupL1)
    #Making sure that dfBestOutGroupOverall is ordered by rownum, this will be important for the next step
    dfBestOutGroupL1 <<- dfBestOutGroupL1[order(dfBestOutGroupL1$rownum),]
    #L2
    dfBestOutGroupL2 <<- by(dfBestOutGroupL2, dfBestOutGroupL2["ind"], head, n=1)
    dfBestOutGroupL2 <<- Reduce(rbind, dfBestOutGroupL2)
    dfBestOutGroupL2 <<- dfBestOutGroupL2[order(dfBestOutGroupL2$rownum),]

    #this will give the correct index value for the rownum column, doing this for both L1 and L2
    j=0:(nrow(dfBestOutGroupL1)-1)
    dfBestOutGroupL1$indexNo <<- (dfBestOutGroupL1$rownum / (nrow(dfAllSeq))-j) * (nrow(dfAllSeq))
    #for merging, need to change ind column in dfAllSeq to numeric first
    indNum <<- with(dfAllSeq, as.numeric(as.character(ind)))
    dfAllSeq$ind <<- indNum
    #using data.tables and setting keys for the merge
    dfAllSeq <<- data.table(dfAllSeq)
    dfBestOutGroupL1 <<- data.table(dfBestOutGroupL1)
    setkey(dfBestOutGroupL1,indexNo)
    setkey(dfAllSeq,ind)
    dfBestOutGroupL1 <<- merge(dfBestOutGroupL1, dfAllSeq, by.x = "indexNo", by.y = "ind", all.x = TRUE)
    #Converting back to dataframe for further use
    dfBestOutGroupL1 <<- as.data.frame(dfBestOutGroupL1)
    dfBestOutGroupL1 <<- dfBestOutGroupL1[order(match(dfBestOutGroupL1[,3],dfMatchOverallLineage1[,26])),]

    #L2
    j=0:(nrow(dfBestOutGroupL2)-1)
    dfBestOutGroupL2$indexNo <<- (dfBestOutGroupL2$rownum / (nrow(dfAllSeq))-j) * (nrow(dfAllSeq))
    dfBestOutGroupL2 <<- data.table(dfBestOutGroupL2)
    setkey(dfBestOutGroupL2,indexNo)
    dfBestOutGroupL2 <<- merge(dfBestOutGroupL2, dfAllSeq, by.x = "indexNo", by.y = "ind", all.x = TRUE)
    dfBestOutGroupL2 <<- as.data.frame(dfBestOutGroupL2)
    dfBestOutGroupL2 <<- dfBestOutGroupL2[order(match(dfBestOutGroupL2[,3],dfMatchOverallLineage2[,26])),]

    #Can rename certain columns to more closely resemble dfMatchOverallLineage dataframes (except the outGroupDist column which would be unique to each lineage of course)
    dfBestOutGroupL1 <<- dfBestOutGroupL1[,c("bin_uri","record_id","values","medianLat","latMin","latMax","binSize","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","medianLon","indexNo","ind")]
    colnames(dfBestOutGroupL1)[3] <<- "outGroupDist"
    dfBestOutGroupL2 <<- dfBestOutGroupL2[,c("bin_uri","record_id","values","medianLat","latMin","latMax","binSize","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","medianLon","indexNo","ind")]
    colnames(dfBestOutGroupL2)[3] <<- "outGroupDist"

    #adding an ingroup pairing column to each bestoutgroup dataframes
    dfBestOutGroupL1$inGroupPairing <<- dfMatchOverallLineage1$inGroupPairing
    dfBestOutGroupL2$inGroupPairing <<- dfMatchOverallLineage2$inGroupPairing

    #Can now merge our outgroups to our associated lineages
    #Each outgroup will be duplicated for each member of the pairing however the distances should be unique to represent each distance from each outgroup to each lineage

    #.x beside headings are the ingroup pairing data (except unique columns), scroll far enough to the right and headings with .y are the outgroup related columns
    dfMatchOverallLineage1 <<- merge(dfMatchOverallLineage1, dfBestOutGroupL1, by.x = "inGroupPairing", by.y = "inGroupPairing")
    dfMatchOverallLineage2 <<- merge(dfMatchOverallLineage2, dfBestOutGroupL2, by.x = "inGroupPairing", by.y = "inGroupPairing")

    #Then we can revise dfMatchOverallBest to correctly reflect both pairings with correct outgroup distances
    #We do this by using rbind to combine both lineages together into the dfMatchOverallBest dataframe
    dfMatchOverallBest <<-  rbind(dfMatchOverallLineage1, dfMatchOverallLineage2)
    #Then order by inGroupPairing once again for a better organization of the dataframe
    dfMatchOverallBest <<- dfMatchOverallBest[order(dfMatchOverallBest$inGroupPairing),]

    #As a last step for outgroup determination, if a suitable outgroup could not be found and is NA for a pairing, then we can filter these pairings out
    noOutGroup<<-which(is.na(dfMatchOverallBest$bin_uri.y))
    if(length(noOutGroup) >0){
      dfMatchOverallBest<<-dfMatchOverallBest[-noOutGroup,]}

    #Also renumber pairings starting from 1 again
    dfMatchOverallBest$inGroupPairing <<- rep(1:(nrow(dfMatchOverallBest)/2), each = 2)
    dfMatchOverallLineage1$inGroupPairing <<- 1:nrow(dfMatchOverallLineage1)
    dfMatchOverallLineage2$inGroupPairing <<- 1:nrow(dfMatchOverallLineage2)
    dfMatchOverallBest
  })
  
  # Create the map
  output$worldmap <- renderLeaflet({
    leaflet() %>%
      addTiles() %>% # Add default OpenStreetMap map tiles
      setView(lng = -93.85, lat = 37.45, zoom = 4)
  })
  
  leafletProxy("worldmap", data = zipdata) %>%
    clearShapes() %>%
    addCircles(~longitude, ~latitude, radius=radius, layerId=~zipcode,
               stroke=FALSE, fillOpacity=0.4, fillColor=pal(colorData)) %>%
    addLegend("bottomleft", pal=pal, values=colorData, title=colorBy,
              layerId="colorLegend")
  })
  
  # Show a popup at the given location (not done yet)
  showLocationPopup <- function(lat, lng) {
    content <- as.character(tagList(
      tags$h4("Name of Phylum:", as.character(colnames(dfInitial)[4])),
      tags$strong(HTML(sprintf("%s, %s",
                               selectedZip$city.x, selectedZip$state.x
      ))), tags$br(),
      sprintf("Longitude: %s", dollar(selectedZip$income * 1000)), tags$br(),
      sprintf("Latitude: %s%%", as.integer(selectedZip$college)), tags$br()
    ))
    leafletProxy("worldmap") %>% addPopups(lng, lat, content)
  }
  
  # When map is clicked, show a popup with city info
  observe({
    leafletProxy("worldmap") %>% clearPopups()
    event <- input$map_shape_click
    if (is.null(event))
      return()
    
    isolate({
      showZipcodePopup(event$id, event$lat, event$lng)
    })
  })
  
  output$url <- DT::renderDataTable( 
    outgroupInput(),
    options = list(scrollX = TRUE)
  )
  
  # output$distPlot <- renderPlot({
  #
  # })


})
