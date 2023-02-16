HUM_STAR=star-index-dir
HUM_BOWTIE_DIR=human-bowtie
HUM_BOWTIE=$(HUM_BOWTIE_DIR)/bowtie-index
KRONA=krona


# other files krona makes, below with $KRONA prefix
# KTAX=division.dmp  gencode.dmp  nodes.dmp       citations.dmp     merged.dmp   readme.txt      taxonomy.tab delnodes.dmp      gc.prt  names.dmp # KACC=all.accession2taxid.sorted
# accession2taxid/ -- gz files in there
# taxdump.tar.gz  

.PHONY: all human krona smallmem

smallmem:  krona $(HUM_BOWTIE).1.bt21

all: human krona

human:  $(HUM_STAR)/SAindex  $(HUM_BOWTIE).1.bt21

krona: $(KRONA)/all.accession2taxid.sorted $(KRONA)/nodes.dmp

$(KRONA)/all.accession2taxid.sorted:
	ktUpdateTaxonomy.sh  --accessions $(KRONA)

$(KRONA)/nodes.dmp: 
	ktUpdateTaxonomy.sh $(KRONA) 

$(HUM_STAR)/SAindex: hg38.fa panTro4.fa
	mkdir -p $(HUM_STAR)
	STAR --runMode genomeGenerate --runThreadN 24 --genomeFastaFiles hg38.fa panTro4.fa --genomeDir $(HUM_STAR) --limitGenomeGenerateRAM 35085737173  2>&1 > star-build.log 

$(HUM_BOWTIE).1.bt21: hg38.fa panTro4.fa
	mkdir -p $(HUM_BOWTIE_DIR)
	bowtie2-build  panTro4.fa,hg38.fa $(HUM_BOWTIE) --large-index 2>&1 > bowtie-build.log 

panTro4.fa: panTro4.fa.gz
	gunzip panTro4.fa.gz

hg38.fa: hg38.fa.gz
	gunzip hg38.fa.gz

hg38.fa.gz:
	wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -N

panTro4.fa.gz:
	wget http://hgdownload.cse.ucsc.edu/goldenPath/panTro4/bigZips/panTro4.fa.gz -N

# TODO: ktUpdateTaxonomy.sh $SOMEDIR
