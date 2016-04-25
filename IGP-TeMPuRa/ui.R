# ui.R

#Packaged required
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

shinyUI(navbarPage("TeMPuЯa", theme = "bootstrap.css", id="nav", position = c("fixed-top"),
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
        div(DT::dataTableOutput("url", width = 1000)),
        #div(style='height:600px; width:950px; overflow:scroll',
        #    DT::dataTableOutput("url", width = 850)),
        br(),
        plotOutput("distPlot")
      )
    )
  ),
  tabPanel("Genetic Distance Models Info",
    h1("Genetic distance models:"),
    a("Link to more explanation for the distance models used in R", href = "http://svitsrv25.epfl.ch/R-doc/library/ape/html/dist.dna.html", target='_blank'),
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
      input from the user in order to perform the comparisons. To look at the source
        code, please go to our Github page below."),
      a('Github', href='https://github.com/IGP-TeMPuRa/web-module', target='_blank'),
      h1("Purpose"),
      p("To develop a method to organize populations by tropical and temperate environments
        in order to analyze the effect of differential rate of evolution by nucleotide
        substitution due to geographical distances for indicators in the environment or
        disease transmitters."),
      h1("Scope"),
      h2("Description of the Data"),
      p("The data for the Project consists of species information (genetic sequence, 
location of capture, etc.), which are currently stored within the Barcode of Life 
Database (BOLD). Historical data is also found within BOLD. To be processed, the 
data files need to be formatted as: .tsv format."),
      h2("Description of Target Users"),
      p("It is expected that the user, who will primarily be processing the data, 
will be a researcher within the field of biodiversity.  The user will not necessarily 
need to have experience with R, but knowledge in manipulating R scripts would be helpful."),
      br(),
      p("To manage the TeMPuЯa system, a researcher familiar with R will need to 
        be available occasionally for general maintenance and post-modification of 
        the system as needs arise."),
      h2("Description of the Proposed Method for Data Processing"),
      p("The proposed process flow for the TeMPuЯa project is as follows:"),
      tags$ul(
        tags$li("User logs onto a computer having access to the internet"),
        tags$li("User chooses the taxonomic group at which to conduct the search"),
        tags$li("The user will need to specify a latitudinal difference. 
                This will tell the program how far or how near to search"),
        tags$li("The user will then specify a genetic similarity threshold.
                This will tell the program to include or exclude pairs of 
                species based on how similar they are according to a 
                predetermined genetic barcode"),
        tags$li("The program will then display all species pairs 
                that fit according to the user-designated specifications on a world map"),
        tags$li("The program will also output relevant information in a .csv file"),
        tags$li("The program will perform relevant statistical tests on pairings"),
        tags$li("The program will plot the pairings as a function of relative outgroup distance")
      ),
      h2("Description of the R Script"),
      p("The R script will consist of 1 file that performs the following functions:"),
      tags$ul(
        tags$li("The file will scrape BOLD for the user designated species and taxonomic level"),
        tags$li("The file will download a .tsv file containing the relevant species 
                data and parse it"),
        tags$li("The file will generate sister pairs of related species and organize 
                them in multiple matrices within R"),
        tags$li("The file will perform analysis on the sister pairs based on the 
                user-specified latitude difference and genetic similarity threshold"),
        tags$li("The file will then perform further statistical tests on the generated pairings")
      ),
      p("The R script itself does not generate a user interface. To facilitate this function, 
        the R script will be hosted as a web application on Shiny."),
      h2("Description of User Interface"),
      p("Graphical user interface (GUI) for R script:"),
      tags$ul(
        tags$li("Taxonomic search: tells the script where to look within BOLD and 
                begin downloading the file"),
        tags$li("Latitudinal difference: slider bar (left side = low latitudinal difference, 
                right side = high latitudinal difference)"),
        tags$li("Genetic similarity threshold: slider bar (left side = low genetic similarity, 
                right side = high genetic similarity)"),
        tags$li("Pairs of species that fit the user's specifications displayed on a world 
                map (colour coded)"),
        tags$li("Table of the results is shown below the world map that can be downloaded"),
        tags$li("Plots for the statistical tests are shown as well")
      )
    ),
    tabPanel("About Us",
      h1("Team"),
      h3("Winfield Ly, Matthew Orton, and David Lee"),
      p("We are a group of budding bioinformaticians/programmers who are classmates in the 
        Bioinformatics Graduate Certificate Program at Seneca College in 2015-2016 to learn 
        and develop our programming skills in Perl, Java, and R. We are fast learners, and 
        work well in a team with excellent communication skills."),
      h1("Special Thanks"),
      p("We would like to thank:"),
      tags$ul(
        tags$li("Our principal investigator, Dr. Sarah Adamowicz,
                 for giving us this opportunity to collaborate between schools 
                 during the school year."),
        tags$li("Keshav Dial and Bilal Athar during the initial stages of this project."),
        tags$li("Monica Wong with help in R"),
        tags$li("Dan Fieldhouse, Program Coordinator of the Bioinformatics Graduate Certificate 
                 Program at Seneca College.")
        )
    )
  )
))
