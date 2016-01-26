#Installation of packages for web scraping in R:
#instructions can be found:http://blog.rstudio.org/2014/11/24/rvest-easy-web-scraping-with-r/
#CRAN page: https://cran.r-project.org/web/packages/rvest/index.html

install.packages("xml2")
install.packages("rvest")
install.packages("httr")
install.packages("selectr")
install.packages("magrittr")
library(rvest)

#Still have to get sequence quality to work properly
#Scrape BIN specimen record (just chose a random one to illustrate):
BINspecimenrecord <- read_html("http://www.boldsystems.org/index.php/Public_RecordView?processid=BLPDY748-11")
BINspecimenrecord %>%
html_node(".subtitles") %>%
html_text()
#[1]"BLPDY748-11"

#Scrape Species name (might come in handy to have the species names stored as well)
SpeciesName <- read_html("http://www.boldsystems.org/index.php/Public_RecordView?processid=BLPDY748-11")
SpeciesName %>%
html_node("tr~ tr+ tr em") %>%
html_text()
#[1] "elachJanzen01 Janzen864"

#Scraping latitude using a random BIN specimen from database
#After using SelectorGadget chrome extension to find html node (most likely we could just use these same nodes for every BIN):
latitude <- read_html("http://www.boldsystems.org/index.php/Public_RecordView?processid=BLPDY748-11")
latitude %>%
html_node(".binDataTable+ .binDataTable tr:nth-child(6) td:nth-child(3)") %>%
html_text() %>%
as.numeric()
#[1] 10.765

#Sequence data:
DNAseq <- read_html("http://www.boldsystems.org/index.php/Public_RecordView?processid=BLPDY748-11")
DNAseq %>%
html_node("tr:nth-child(5) pre") %>%
html_text() 
#[1] "AACATTATATTTTATTTTTGGAATTTGAGCAGGTATAGTCGGAACTTCTTTAAGATTATTAATTCGAGCTGAATT\nAGGTAACCCAGGGTCTTTAATTGGAGATGATCAAATTTATAATACAATTGTTACA
#GCTCATGCTTTTATTATAAT\nTTTTTTTATAGTAATACCTATTATAATTGGAGGATTTGGAAATTGACTTGTTCCTTTAATATTAGGAGCTCCTGA\nTATAGCTTTCCCCCGAATAAATAACATAAGTTTCTGAT
#TATTGCCCCCTTCTCTTACTCTATTAATTTCTAGAAG\nAATTGTAGAGAATGGGGCAGGGACTGGATGAACAGTTTACCCCCCACTTTCATCTAATATTGCCCATGGGGGTAG\nATCAGTAGATTTAGCAATTTT
#TTCTTTACATTTAGCTGGTATTTCTTCAATTTTAGGGGCAATCAATTTTATTAC\nTACTATTATTAATATACGATTAAATAATATATCATTTGATCAATTACCTTTATTTGTTTGAGCAGTGGGAATTAC\nAGCT
#TTATTATTACTTTTATCTTTACCTGTTTTAGCTGGAGCTATCACTATATTATTAACTGATCGAAATTTAAA\nTACATCTTTTTTTGATCCTGCTGGAGGAGGAGATCCTATTTTATACCAACATTTATTT"
