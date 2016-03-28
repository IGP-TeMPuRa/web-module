# server.R

library(shiny)
library(leaflet)
library(foreach)
library(ape)
library(data.table)
library(DT)
source("tsvtoDataFrame.R")

shinyServer(function(input, output, session) {
  
  observe({
    query <- parseQueryString(session$clientData$url_search)
    if (!is.null(query[['taxonomy']])) {
      updateTextInput(session, "taxonomy", value = query[['taxonomy']])
    }
  })
  
  observe({
    query <- parseQueryString(session$clientData$url_search)
    if (!is.null(query[['geography']])) {
      updateTextInput(session, "geography", value = query[['geography']])
    }
  })
  
  observe({
    query <- parseQueryString(session$clientData$url_search)
    if (!is.null(query[['latitude']])) {
      updateSliderInput(session, "latitude", value = query[['latitude']])
    }
  })
  
  observe({
    query <- parseQueryString(session$clientData$url_search)
    if (!is.null(query[['genetic']])) {
      updateSliderInput(session, "genetic", value = query[['genetic']])
    }
  })
  
  observe({
    query <- parseQueryString(session$clientData$url_search)
    if (!is.null(query[['outgroups']])) {
      updateSliderInput(session, "outgroups", value = query[['outgroups']])
    }
  })
  
  observe({
    query <- parseQueryString(session$clientData$url_search)
    if (!is.null(query[['distanceModels']])) {
      updateSliderInput(session, "distanceModels", value = query[['distanceModels']])
    }
  })
  
  # Create the map
  output$worldmap <- renderLeaflet({
    leaflet() %>%
      addTiles() %>% # Add default OpenStreetMap map tiles
      setView(lng = -93.85, lat = 37.45, zoom = 4) 
  })
  
  # Store the user input into separate variables
  #taxonInput <- input$taxonomy
  #geographyInput <- input$geography
  #latitudeInput <- input$latitude
  #geneticInput <- input$genetic
  #outgroupInput <- input$outgroups
  #distanceInput <- input$distanceModels
  
  
  
  textInput <- reactive({
    var1 <- "http://www.boldsystems.org/index.php/API_Public/combined?taxon="
    var2 <- "&geo="
    var3 <- "&format=tsv"
    paste(c(var1), c(input$taxonomy), c(var2), c(input$geography), c(var3), collapse = '')
    
  })
  
  output$text <- renderText({
    textInput
  })
  
  #observe({
  #  leafletProxy("worldmap", data = filteredData) %>%
  #    clearShapes() %>%
  #    addMarkers("worldmap", lng = NULL, lat = NULL, popup = NULL, options = markerOptions,
  #               data = getMapData())
  #})
  
  output$url <- DT::renderDataTable( 
    dfMatchOverallBest,
    options = list(scrollX = TRUE)
  )


})