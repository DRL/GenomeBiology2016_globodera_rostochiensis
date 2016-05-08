# GenomeBiology2016_globodera_rostochiensis
Scripts used in "The genome of the yellow potato cyst nematode, Globoderia rostochiensis, reveals insights into the bases of parasitism and virulence"

- extractRegionFromCoordinates.py
```
Usage: ./extractRegionFromCoordinates.py [GFF] [FASTA] [US] [DS] [UE] [DE] 

	[GFF] : Intron features have to be present in GFF (use Genometools)
	[US] : Positions upstream of start of intron feature in GFF
	[DS] : Positions downstream of start of intron feature in GFF
	[UE] : Positions upstream of end of intron feature in GFF
	[DS] : Positions downstream of end of intron feature in GFF

   - Extracting splice sites : 

   ./extractRegionFromCoordinates.py nGr.v1.0.gff3 nGr.v1.0.fa 0 1 1 0 
```
- rbbh.py
```
Usage: ./rbbh.py [BLAST1] [FASTA1] [BLAST2] [FASTA2] [EVAL] [COV]

			[BLAST1] : BLAST result of G. pallida vs. G. rostochiensis
			[FASTA1] : Protein sequences of G. pallida
			[BLAST2] : BLAST result of G. rostochiensis vs. G. pallida
			[FASTA2] : Protein sequences of G. rostochiensis
			[EVAL]   : Float of minimal evalue theshold of BLAST hit
			[COV]    : Float of minimal coverage threshold of BLAST hit
```
- parse_iadhore.py
```
Usage: ./parse_iadhore.py [GROS_BED] [GROS_LENGTH] [GPAL_BED] [GPAL_LENGTH] [ANCHORPOINTS] [BASECLUSTERS] [GROS_ASSEMBLY] [GPAL_ASSEMBLY] 

	[GROS_BED] 	: Bed file of G. rostochiensis
	[GROS_LENGTH] 	: Length file of G. rostochiensis
	[GPAL_BED] 	: Bed file of G. pallida
	[GPAL_LENGTH] 	: Length file of G. pallida
	[ANCHORPOINTS] 	: Anchorpoints file provided by i-ADHoRe
	[BASECLUSTERS] 	: Baseclusters file provided by i-ADHoRe
	[GROS_ASSEMBLY] : Assembly file of G. rostochiensis
	[GPAL_ASSEMBLY] : Assembly file of G. pallida
```
- parse_snpeff.py
```
Usage : ./parse_snpeff.py [VCF.gz]

            [VCF.gz]  : Tabix indexed SNPeff VCF file
```
