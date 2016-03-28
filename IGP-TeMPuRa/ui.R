# ui.R

library(shiny)
library(leaflet)
library(foreach)
library(ape)
library(data.table)
library(DT)

# Choices for the genetic distance model
geneticDistanceModel <- c(
  "raw" = "raw", 
  "JC69" = "JC69", 
  "K80" = "K80", 
  "F81" = "F81", 
  "K81" = "K81", 
  "F84" = "F84", 
  "BH87" = "BH87",
  "T92" = "T92",
  "TN93" = "TN93",
  "GG95" = "GG95",
  "logdet" = "logdet",
  "paralin" = "paralin"
)

shinyUI(navbarPage("TeMPuЯa", id="nav", position = c("fixed-top"),
  # needed to keep fixed-top navbar from obscuring content
  header = tags$style(type = "text/css", "body {padding-top: 70px;}"),
  collapsible = "true",
  tabPanel("Tool",
    h1("Instructions"),
    p("Placeholder"),
    sidebarLayout(
      sidebarPanel(
        textInput("taxonomy", label = h4("Enter taxonomy group:"), value = "Porifera"),
        textInput("geography", label = h4("Enter geographical location:"), value = "all"),
        sliderInput("latitude", label = h4("Latitude difference"), min = 10, max = 30, value = 20),
        sliderInput("genetic", label = h4("Genetic similarity threshold"), min = 10, max = 20, value = 15),
        sliderInput("outgroups", label = h4("Select a distance from the outgroup"), min = 1, max = 2, value = 1.3, step = 0.1),
        selectInput("distanceModels", label = h4("Select a genetic distance model"), geneticDistanceModel, selected = "K80"),
        submitButton("Submit"),
        br(),
        downloadButton("download", label = "Download CSV")
      ),
      mainPanel(
        leafletOutput("worldmap"),
        br(),
        div(style='height:300px; width:850px; overflow:scroll',
            DT::dataTableOutput("url", width = 850)),
        textOutput("text")
      )
    )
  ),
  tabPanel("Genetic Distance Models Info",
    h1("Genetic distance models:"),
    a("Link to more explanation for the distance models used in R", href = "http://svitsrv25.epfl.ch/R-doc/library/ape/html/dist.dna.html"),
    br(),
    p(strong("raw:") ,"This is simply the proportion or the number of sites that differ between each pair of sequences. This may be useful to draw 'saturation plots'."),
    p(strong("JC69:") ,"This model was developed by Jukes and Cantor (1969)."),
    p(strong("K80:") ,"The distance derived by Kimura (1980), sometimes referred to as 'Kimura's 2-parameters distance'."),
    p(strong("F81:") ,"Felsenstein (1981) generalized the Jukes-Cantor model."),
    p(strong("K81:") ,"This model is called the Kimura's 'three substitution types model' (3ST), and is sometimes referred to as 'Kimura's 3-parameters distance'."),
    p(strong("F84:") ,"This model generalized K80, and was first introduced by Felsenstein in 1984."),
    p(strong("BH87:") ,"Barry and Hartigan (1987)."),
    p(strong("T92:") ,"Tamura (1992) generalized the Kimura model."),
    p(strong("TN93:") ,"Tamura and Nei (1993) model."),
    p(strong("GG95:") ,"Galtier and Gouy (1995) model."),
    p(strong("logdet:") ,"The Log-Det distance, developed by Lockhart et al. (1994), is related to BH87. However, this distance is symmetric."),
    p(strong("paralin:") ,"Lake (1994) developed the paralinear distance which can be viewed as another variant of the Barry-Hartigan distance.")
  ),
  navbarMenu(
    "About",
    tabPanel("More Information",
      h1("Description of this Tool"),
      p("TeMPuЯa is an R pipeline that can effectively perform phylogenetic comparisons on 
      large numbers of species based on latitude and genetic similarity. This page accepts
      input from the user in order to perform the comparisons.")       
    ),
    tabPanel("About Us",
      h1("Team"),
      h3("Winfield Ly"),
      p("Stuff about you"),
      h3("Matthew Orton"),
      p("Stuff about you"),
      h3("David Lee"),
      p("Stuff about you")
    )
  )
))
