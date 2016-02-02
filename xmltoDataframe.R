#Must use XMl package for parsing
require("XML")

#Should work at any taxonomic level using the BOLD API!
#Just using brachipoda phylum as an example since it has a smaller number of members to it
#make sure you type format as xml, taken from the BOLD API webpage: http://www.boldsystems.org/index.php/resources/api?type=webservices
#also only specifying the COI-5P sequence for now

#first parsing the API xml file
#Could be any taxonomic level where it says "taxon=" (phylum, order, family etc.) and any country or continent where it says "geo=". For global, simply type "all" instead of a specific region.
xml1 = xmlParse("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Brachiopoda&geo=all&marker=COI-5P&format=xml")

#putting everything into a large dataframe for now
df = xmlToDataFrame(xml1)

#Comparison package will determine which record id's have lat and lon values since we dont need ids that dont have this info
install.packages("compare")
library(compare)

#Latitude and longitude dataframes
dflon <- xmlToDataFrame(getNodeSet(xml1, "//lon"))
dflat <- xmlToDataFrame(getNodeSet(xml1, "//lat"))

#Sequence datafame
seq <- xmlToDataFrame(getNodeSet(xml1, "//nucleotides"))

# Load the package required to read XML files.
library("XML")
