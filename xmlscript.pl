use strict;
use warnings;
use XML::LibXML;

open my $xml, '<', 'bold_data2.xml';
binmode $xml;

my $file = XML::LibXML->load_xml(IO => $xml);
#my $file2 = XML::LibXML->load_xml(IO => $xml);
$file =~ s/\n\s*//g;
my @array;
my $name;
# $name = $1 if
while ($file =~ /(?<=bin_uri>BOLD:)(.*?)(?=<\/bin_uri)/g){
    push @array, $1;
}

foreach (@array){
    print "$_\n";
}
    
