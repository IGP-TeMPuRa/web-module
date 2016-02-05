# THIS IS A VERSION OF PARSING DATA AND PREPARING THE DATA.FRAME FOR R USING PERL REGEX /// ALSO CONTAINS QUALITY ASSESSMENT OF INFORMATION
# 
#
# MODULE METHOD APPEARS TO BE SIMPLER /// CAN STILL USE THIS FOR TESTING AGAINST MODULE METHOD OF CHOICE.
#
# VALUE
# Bilal Athar

use strict;
use warnings;
use XML::LibXML;

my (@binArray, @binContent, $counter, $temp, @speciesName, @nucleotides, @lattitude, @longitude, $tempRegion, $holdlatlon, $lengthSeq, $percentWild, $counter2, $seqHolder, $substrOutput, $errorCount);
$errorCount = 0;

open my $xml, '<', 'bold_data(2).xml';
binmode $xml;
my $file = XML::LibXML->load_xml(IO => $xml);

while ($file =~ /(?<=bin_uri>)(.*?)(?=<\/bin_uri)/g) {
    push(@binArray, $1);    
}

while ($file =~ /(?<=\/bin_uri>)(.*?)(?=<\/record)/gs) {
    push(@binContent, $1);
    $temp = $1;
    if ($temp =~ /(<species(.*?)name>(.*?)<\/name)/gs) {
        push(@speciesName, $3);
    } 
    else {
        push(@speciesName, "");
    }
    if ($temp =~ /(?<=nucleotides>)(.*?)(?=<\/nucleotides)/g) {
        push(@nucleotides, $1);
    }
    else {
        push(@nucleotides, "");
    }
    if ($temp =~ /(?<=<lon>)(.*?)(?=<\/lon>)/gs) {
        push(@longitude, $1);
        if ($temp =~ /(?<=<lat>)(.*?)(?=<\/lat>)/gs) {
            push(@lattitude, $1);
        }
        else {
            push(@lattitude, "");
        }
    }
    elsif ($temp =~ /(?<=<region>)(.*?)(\d.*)(?=<\/region>)/g) {
        $tempRegion = "$2";
        if ($tempRegion =~ /((\d.*\s)([N|S]))/) {
            if ($3 eq "N") {
                push (@lattitude, $2);
            }
            else {
                $holdlatlon = "-";
                $holdlatlon .= "$2";
                push (@lattitude, $holdlatlon);
            }         
        }
        if ($tempRegion =~ /(?<=,\s)(\d.*\s)(E|W)/) {
            if ($2 eq "E") {
                push (@longitude, $1);                
            }
            else {
                $holdlatlon = "-";
                $holdlatlon .= "$1";
                push (@longitude, $holdlatlon);
            }
        }
    } else {
        push (@lattitude, "");
        push (@longitude, "");
    } 
}

for ($counter=0; $counter<scalar(@binArray); $counter++) {
    if ($binArray[$counter] eq "" || $speciesName[$counter] eq "" || $nucleotides[$counter] eq "" || $lattitude[$counter] eq "" || $longitude[$counter] eq "") {
        splice (@binArray, $counter, 1);
        splice (@speciesName, $counter, 1);
        splice (@nucleotides, $counter, 1);
        splice (@lattitude, $counter, 1);
        splice (@longitude, $counter, 1);
        $counter--;    
    }   
}

for ($counter=0; $counter<scalar(@binArray); $counter++) {
    $seqHolder = $nucleotides[$counter];
    $lengthSeq = length($seqHolder);  
    if ($lengthSeq < 500) {
        splice (@binArray, $counter, 1);
        splice (@speciesName, $counter, 1);
        splice (@nucleotides, $counter, 1);
        splice (@lattitude, $counter, 1);
        splice (@longitude, $counter, 1);
        $counter--;
    }
    else {               
        for ($counter2 = 0; $counter2 < $lengthSeq; $counter2++) {
            $substrOutput = substr($seqHolder, $counter2, 1);
            if ($substrOutput eq "N") {
                $errorCount++;
            }        
        }
        $percentWild = ($errorCount / $lengthSeq);
        if ($percentWild > 0.01) {
            splice (@binArray, $counter, 1);
            splice (@speciesName, $counter, 1);
            splice (@nucleotides, $counter, 1);
            splice (@lattitude, $counter, 1);
            splice (@longitude, $counter, 1);
            $counter--;
        }
    }    
}

my $counter3 = 0;
print "\n\n\nBINURI's\n\n\n";

foreach (@binArray) {
    print "$counter3:  $_\n";
    $counter3++;
}

print "\n\n\nSPECIES\n\n\n";

$counter3 = 0;
foreach (@speciesName) {
    print "$counter3:  $_\n";
    $counter3++;
}

print "\n\n\nLATITTUDE\n\n\n";

$counter3 = 0;
foreach (@lattitude) {
    print "$counter3:  $_\n";
    $counter3++;
}

print "\n\n\nLONGITUDE\n\n\n";

$counter3 = 0;
foreach (@longitude) {
    print "$counter3:  $_\n";
    $counter3++;
}

print "\n\n\nSEQUENCE\n\n\n";

$counter3 = 0;
foreach (@nucleotides) {
    print "$counter3:  $_\n";
    $counter3++;
}

for ($counter=0; $counter<scalar(@binArray); $counter++) {
    print "BIN URI: ", $binArray[$counter], " "x5;
    print "SPECIES NAME: ", $speciesName[$counter], " "x5;
    print "LONGITUDE: ", $longitude[$counter], " "x5;
    print "LATTITUDE: ", $lattitude[$counter], " "x5;
    # print "SEQUENCE: ", $nucleotides[$counter], " "x5;
    print "\n";    
}

$counter3--;
print "\n\n\tThis DataSet has $counter3 rows of Complete and Satisfactory Data\n\n";




   
#- EXTRACTION -#
    #- BIN UI's into an array ALL of them sequentially.
    #- Array of all Bin UI supplementary information
    #- Array of all matched extracted supplementary info from Array 2. ( Longitude, Lattitude, Sequence, Species Name)


#- SCREENING FOR QUALITY -#
    #- Sequence Missing or Not, Longitude and Latitude missing or NOT, % of Sequence thats N or -, 
    #- Conversion of Array information into seperate Arrays Ready for R Input data.frame list format
    #- Removal of invaluable data point index from each Array.


#- QUESTIONS -#
    #- WHAT TO DO WHEN NO LON/LAT BUT INSTEAD A COUNTRY -#
