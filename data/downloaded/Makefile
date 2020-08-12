all: timestamp.txt

timestamp.txt: epifactors.genes.csv.gz epifactors.complexes.csv.gz  biogrid.homo.sapiens.tab2.txt.gz mygene/cache.h5cache ebi.complexes.tsv.gz
	date > timestamp.txt

epifactors.genes.csv.gz:
	wget --user-agent="Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36" http://epifactors.autosome.ru/public_data/current/EpiGenes_main.csv -O epifactors.genes.csv
	gzip epifactors.genes.csv

epifactors.complexes.csv.gz:
	wget --user-agent="Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36" http://epifactors.autosome.ru/public_data/current/EpiGenes_complexes.csv -O epifactors.complexes.csv
	gzip epifactors.complexes.csv

ebi.complexes.tsv.gz:
	wget ftp://ftp.ebi.ac.uk/pub/databases/intact/complex/2019-07-19/complextab/homo_sapiens.tsv -O ebi.complexes.tsv
	gzip ebi.complexes.tsv

biogrid.homo.sapiens.tab2.txt.gz:
	wget https://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.5.174/BIOGRID-ORGANISM-3.5.174.tab2.zip -O biogrid.organism.tab2.zip
	unzip -u biogrid.organism.tab2.zip -d biogrid.organism
	gzip biogrid.organism/BIOGRID-ORGANISM-Homo_sapiens-3.5.174.tab2.txt
	mv biogrid.organism/BIOGRID-ORGANISM-Homo_sapiens-3.5.174.tab2.txt.gz biogrid.homo.sapiens.tab2.txt.gz
	rm -rf biogrid.organism/
	rm biogrid.organism.tab2.zip

mygene/cache.h5cache: mygene/mygene.list
	python mygene/populate_cache.py

.PHONY: all
