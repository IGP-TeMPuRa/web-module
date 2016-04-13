# ui.R

library(shiny)
library(leaflet)
library(foreach)
library(ape)
library(data.table)
library(DT)
library(nlme)

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
    p("The left sidebar panel contains the inputs that the user can enter into the pipeline
      to fit for their analysis. The taxonomy group must already exist on the BOLD database. On the right side of the side panel is the map that will 
      show the markers of where the species is located in the world. Below the map will show
      the results in table format and the statistical analysis of the data. The top of the page
      links to more pages that has more information about this project, about the genetic 
      distance models used in R, and about us."),
    sidebarLayout(
      sidebarPanel(
        textInput("taxonomy", label = h4("Enter taxonomy group:"), value = "Porifera"),
        textInput("geography", label = h4("Enter geographical location:"), value = "all"),
        sliderInput("latitude", label = h4("Latitude difference"), min = 10, max = 30, value = 20),
        sliderInput("genetic", label = h4("Genetic similarity threshold"), min = 0.10, max = 0.20, value = 0.15),
        sliderInput("outgroups", label = h4("Select a distance from the outgroup"), min = 1, max = 2, value = 1.3, step = 0.1),
        selectInput("distanceModels", label = h4("Select a genetic distance model"), geneticDistanceModel, selected = "TN93"),
        submitButton("Submit"),
        br(),
        downloadButton("download", label = "Download CSV")
      ),
      mainPanel(
        leafletOutput("worldmap"),
        br(),
        div(style='height:600px; width:850px; overflow:scroll',
            DT::dataTableOutput("url", width = 850)),
        br(),
        plotOutput("distPlot")
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
      large numbers of species based on latitude and genetic similarity. This website accepts
      input from the user in order to perform the comparisons."),
      h1("Purpose"),
      p("To develop a method to organize populations by tropical and temperate environments
        in order to analyze the effect of differential rate of evolution by nucleotide
        substitution due to geographical distances for indicators in the environment or
        disease transmitters."),
      h1("Scope"),
      p("The scope of this project is to design and construct an analytical software tool for
        Dr. Sarah Adamowicz in the Centre for Biodiversity Genomics at the University
        of Guelph. This tool will accept the species and criteria from the user, obtain
        the required information from the Barcode of Life Data Systems (BOLD), screen
        the species according to the quality of the sequence of the genetic marker, GPS
        coordinates and the genetic distance, and output both the pairs of sister lineages
        and outgroup sequences in CSV format suitable for follow-on analysis in the
        software PAML. All of the data are properties of BOLD and the University of Guelph.")
    ),
      h3("Link to the source code"),
      a("Github", href="https://github.com/IGP-TeMPuRa/web-module", target="_blank"),
    tabPanel("About Us",
      h1("Team"),
      h3("Winfield Ly, Matthew Orton, and David Lee"),
      p("We are a group of budding bioinformaticians/programmers who are classmates in the 
        Bioinformatics Graduate Certificate Program at Seneca College in 2015-2016 to learn 
        and develop our programming skills in Perl, Java, and R. We are fast learners, and 
        work well in a team with excellent communication skills."),
      h3("Special Thanks"),
      p("We would like to thank Keshav Dial and Bilal Athar during the initial stages of this 
        project. We would also like to thank our PI, Dr. Sarah Adamowicz, for giving us this
        opportunity to collaborate between schools during the school year. Finally, we would
        like to thank Dan Fieldhouse, Program Coordinator of the Bioinformatics 
        Graduate Certificate Program at Seneca College.")
    )
  )
))
