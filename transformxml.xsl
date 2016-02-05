<xsl:transform xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
<xsl:output version="1.0" encoding="UTF-8" indent="yes" />
<xsl:strip-space elements="*"/>

  <xsl:template match="/">
    <root>
      <xsl:apply-templates select="*"/>
    </root>
  </xsl:template>  

  <xsl:template match="record">
    <xsl:copy>
      <xsl:copy-of select="record_id"/>
      <xsl:copy-of select="bin_uri"/>
	  <xsl:copy-of select="taxonomy/species/taxon/name"/>
      <xsl:copy-of select="collection_event/coordinates/lat"/>
      <xsl:copy-of select="collection_event/coordinates/lon"/>
      <xsl:copy-of select="sequences/sequence/nucleotides"/>
    </xsl:copy>
  </xsl:template>

</xsl:transform>
