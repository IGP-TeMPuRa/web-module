  #This script along with transfomxml.xsl will take the original xml and greatly simplify to the following structure so all nodes are at the same level making operation in perl 
  #and import to R really easy:

  # <xsl:template match="record">
    # <xsl:copy>
      # <xsl:copy-of select="record_id"/>
      # <xsl:copy-of select="bin_uri"/>     
      # <xsl:copy-of select="taxonomy/species/taxon/name"/>
      # <xsl:copy-of select="taxonomy/phylum/taxon/name"/>
      # <xsl:copy-of select="taxonomy/class/taxon/name"/>
      # <xsl:copy-of select="taxonomy/family/taxon/name"/>
      # <xsl:copy-of select="taxonomy/subfamily/taxon/name"/>
      # <xsl:copy-of select="taxonomy/genus/taxon/name"/>
      # <xsl:copy-of select="collection_event/coordinates/lat"/>
      # <xsl:copy-of select="collection_event/coordinates/lon"/>
      # <xsl:copy-of select="sequences/sequence/nucleotides"/>
    # </xsl:copy>
  # </xsl:template>

  #must use CPAN to install these packages
  #perl -MCPAN -e "install ..."
  
  use XML::LibXSLT;
  use XML::LibXML;

  my $parser = XML::LibXML->new();
  my $xslt = XML::LibXSLT->new();

  #xml file, perl file and xsl script must all be in the same directory
  #insert whichever xml you want and perl will parse it
  my $source = $parser->parse_file('C:\Users\Winfield\Desktop\BOLD_database_xml\bold_data.xml');
  
  #xsl script that simplifies the xml file
  my $style_doc = $parser->parse_file('C:\Users\Winfield\Desktop\BOLD_database_xml\transformxml.xsl');
  my $stylesheet = $xslt->parse_stylesheet($style_doc);

  my $results = $stylesheet->transform($source);

  #to export to XML for use in R, will output to directory perl script is located in
	
  open XML, ">revisedBOLD2.xml";
  print XML $stylesheet->output_string($results);
  close XML;
