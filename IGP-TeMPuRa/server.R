# server.R

#Packages required
library(shiny)
library(leaflet)
library(DT)
library(magrittr)
#install.packages("nlme")
library(nlme)
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
source("http://bioconductor.org/biocLite.R")
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
  
  #source('tsvtoDataFrame.R', local=TRUE)
  source('demo.R', local=TRUE)
  
  #Used to get the latitude and longitude points and convert them back to their original value
  points <- eventReactive(input$latitude, {
    dfMatchOverallBest$medianLat.x <<- dfMatchOverallBest$medianLat.x - 90
    dfMatchOverallBest$medianLon.x <<- dfMatchOverallBest$medianLon.x - 180
    cbind(dfMatchOverallBest$medianLat.x, dfMatchOverallBest$medianLon.x)
  }, ignoreNULL = FALSE)
  
  #content <- paste0(sep = "<br/>","<b>Porifera</b>")
  
  # Create the map
  output$worldmap <- renderLeaflet({
    leaflet() %>%
      addTiles() %>% # Add default OpenStreetMap map tiles
      setView(lng = -93.85, lat = 37.45, zoom = 1) %>%
      addMarkers(data = points()) #%>%
      #addPopups(dfMatchOverallBest$medianLat.x, dfMatchOverallBest$medianLon.x, content,
                #options = popupOptions(closeButton = TRUE))
  })
  
  #Returns a table from the R script to render it on the UI side 
  output$url <- DT::renderDataTable( 
    dfMatchOverallBest,
    options = list(scrollX = TRUE)
  )
  
  #Returns the graph from the R script to render it on the UI side 
  output$distPlot <- renderPlot({
    plot
  })
  
  # downloadHandler() takes two arguments, both functions.
  # The content function is passed a filename as an argument, and
  # it should write out data to that filename.
  output$download <- downloadHandler(
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() { paste0(input$taxonomy, '.csv', sep='') },
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file){
      write.csv(dfMatchOverallBest, file)
    }
  )
})
