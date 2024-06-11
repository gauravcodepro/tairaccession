# tairaccession

- tairaccession: This allows you to analyze the tair accession ids easily and you can see the sample datasets for each of the import function with in the tests directory and the corresponding format are available from [TAIR](https://www.arabidopsis.org). This package has added utilities for analyzing also the Phytozome datasets and also you can plot the desried genes of interest. You can also analyze [Phytozome Araport](https://phytozome-next.jgi.doe.gov) also using this package.The package is under release at PyPI package repository. I have updated this package with additional support for the phytozome in addition to the tair. In an update to this package, few functions on plotting the coding regions, gene regions and exons have been added.There are additional functions such ```prepareFunctionalNamePhytozome```and ```preparegeneNamePhytozome```which will automatically prepare the files as per the release of the phytozome and tair.

- The web documentation is located at [tairaccession](https://gauravcodepro.github.io/tairaccession)

Gaurav Sablok, \
Academic Staff Member,Bioinformatics, \
Institute for Biochemistry and Biology, \
University of Potsdam
Potsdam,Germany 

If you have any questions, please contact at gaurav.sablok@uni-potsdam.de.
## Installation

```bash
$ pip install tairaccession
import tairaccession
print(tairaccession.__version__)
```

## Usage

`tairaccession` can be used to access the tairIDs and also used for the conversion of the tair data and obtaining the relative information from the tair accessions aka AGIs specificed in a file. In case of the genes, you need to specify the gene IDs, whereas in case of the exons, it can take the splice variants also. In case you have the splice variant also, and you want the gene coordinates then also it is feasible as it will strip the variant automatically, while searching for the gene. The tair accession has the following options:


```python
from tairaccession.tairaccession import uniprotAgi
uniportAgi(agi_file, ids_file) 
# agi_file: path to the agi_file 
#ids_file: path to the ids_file
uniprotAgi("/Users/gauravsablok/Desktop/release/Uniprot2AGI-Jul2023.txt", \
                      "/Users/gauravsablok/Desktop/release/test_ids.txt")
[('A0A654EWS4', ['AT1G64990.2', 'AT1G64990.1']),
 ('A0A654F9L3', ['AT3G22240.1', 'AT3G22240.1']),
 ('Q9SRQ8', ['At3g03550']),
 ('Q9LRY7', ['At3g24715.1', 'AT3G24715.2', 'AT3G24715.3']),
 ('Q9FT69', ['At5g27680.1', 'AT5G27680.2', 'AT5G27680.3', 'AT5G27680.4'])]                      
help(uniprotAgi) # for detailed information
````
```python
from tairaccession.tairaccession import agiUniprot
agiUniprot(agi_uniprot_file, ids_file)
# agi_uniprot_file: path to the agi_file 
#ids_file: path to the ids_file
agiUniprot("/Users/gauravsablok/Desktop/release/AGI2uniprot-Jul2023.txt", \
                              "/Users/gauravsablok/Desktop/release/test_ids.txt")
[('AT5G27680.3', ['Q9FT69']),
 ('AT1G64990.1', ['A0A654EWS4']),
 ('AT3G22240.1', ['A0A654F9L3']),
 ('AT3G24715.2', ['Q9LRY7'])]                              
help(agiUniprot) # for detailed information
````
```python
from tairaccession.tairaccession import uniprotTair
uniprotTair(tair_uniprot_file, ids_file) 
# tair_uniprot_file: path to the agi_file 
#ids_file: path to the ids_file
uniprotTair("/Users/gauravsablok/Desktop/release/TAIR2UniprotMapping-Jul2023.txt", \
                                    "/Users/gauravsablok/Desktop/release/test_ids.txt")
[('Q9FT69', ['locus:2180255', 'AT5G27680']),
 ('Q9XIP7', ['locus:2010796', 'AT1G64990']),
 ('Q9LHJ3', ['locus:2091623', 'AT3G22240']),
 ('Q9LRY7', ['locus:2093146', 'AT3G24715']),
 ('Q9SRQ8', ['locus:2096444', 'AT3G03550']),
 ('A0A1P8BBY0', ['locus:2180255', 'AT5G27680'])]  
help(uniprotTair) # for detailed information
````
```python
from tairaccession.tairaccession import agiSpliceCoordinates
agiCoordinates(ids_file, gff_file, gene_type = None)  
#ids_file: path to the ids_file #gff_file: path to the gff file 
#gene_type: this is a #keyworded argument to the agi_coordinates. 
#if gene_type is gene: it will return the gene coordinates for the agi in the file, 
#if gene_type #is exon: it will return the exon coordinates for the agi in the file,
#if gene_type is three_prime_UTR: it will return the three_prime_UTR 
#if gene_type is five_prime_UTR: it will return the five_prime_UTR coordinates for the agi in the file,
#if gene_type is cds: it will return the cds for the agi in the file,
agiCoordinates("/Users/gauravsablok/Desktop/release/test_ids.txt", \
                    "/Users/gauravsablok/Desktop/release/TAIR10_GFF3_genes.gff", \
                                                                gene_type = "gene")
# it provides the strand information also which can be directly use in the jbrowse,
# 1 for the positive strand and the -1 for the negative strand.
[['AT1G01010', 3631, 5899, 1],
 ['AT1G01020', 5928, 8737, -1]]

agiCoordinates("/Users/gauravsablok/Desktop/release/test_ids.txt", \
                    "/Users/gauravsablok/Desktop/release/TAIR10_GFF3_genes.gff", \
                                                                  gene_type = "exon")
# it provides the strand information also which can be directly use in the jbrowse,
# 1 for the positive strand and the -1 for the negative strand.
[['AT1G01010.1', 3631, 3913, 1],
 ['AT1G01010.1', 3996, 4276, 1],
 ['AT1G01010.1', 4486, 4605, 1]]

 agiCoordinates("/Users/gauravsablok/Desktop/release/test_ids.txt", \
                     "/Users/gauravsablok/Desktop/release/TAIR10_GFF3_genes.gff", \
                                                       gene_type = "three_prime_UTR")
# it provides the strand information also which can be directly use in the jbrowse,
# 1 for the positive strand and the -1 for the negative strand.
[['AT1G01010.1', 5631, 5899, 1],
 ['AT1G01020.1', 6437, 6914, -1]]
agiCoordinates("/Users/gauravsablok/Desktop/release/test_ids.txt", \
                     "/Users/gauravsablok/Desktop/release/TAIR10_GFF3_genes.gff", \
                                                       gene_type = "five_prime_UTR") 
# it provides the strand information also which can be directly use in the jbrowse,
# 1 for the positive strand and the -1 for the negative strand.
[['AT1G01010.1', 3631, 3759, 1],
 ['AT1G01010.1', 3631, 3759, 1],
 ['AT1G01020.1', 8667, 8737, -1],
 ['AT1G01020.1', 8667, 8737, -1],
 ['AT1G01020.1', 8667, 8737, -1]]
agiSpliceCoordinates("/Users/gauravsablok/Desktop/CodeTest/release/all_ids.txt", \
                    "/Users/gauravsablok/Desktop/CodeTest/release/TAIR10_GFF3_genes.gff", \
                                                                gene_type = "cds")
[['AT1G01010.1', 3760, 3913, 1],
 ['AT1G01010.1', 3760, 3913, 1],
 ['AT1G01010.1', 3996, 4276, 1],
 ['AT1G01010.1', 3996, 4276, 1],
 ['AT1G01010.1', 4486, 4605, 1]]
help(agiCoordinates) # for detailed information
````

```python
from tairaccession.tairaccession import getagiGO
getagiGO(ids_file, association_file) 
#ids_file: path to the ids_file, #association_file: path to the association file
#provides the information on the AGI and the corresponding GO information.
getagiGO("/Users/gauravsablok/Desktop/release/test_ids.txt", \
                         "/Users/gauravsablok/Desktop/release/gene_association.tair")
('AT5G27680', 'GO:0005634'),
 ('AT5G27680', 'GO:0005694'),
 ('AT5G27680', 'GO:0005737'),
 ('AT5G27680', 'GO:0006268'),
 ('AT5G27680', 'GO:0006310'),
 ('AT5G27680', 'GO:0009378'),
 ('AT5G27680', 'GO:0043138'),
 ('AT5G27680', 'GO:0006281'),
 ('AT5G27680', 'GO:0032508'),
 ('AT5G27680', 'GO:0043138')                         
help(getagiGO) # for detailed information
````
```python
from tairaccession.tairaccession import  getagiTAIR
 getagiTAIR(ids_file, association_file) 
 #ids_file: path to the ids_file, 
 #association_file: path to the association file 
 #provides the information on the AGI and the corresponding TAIR communication.
 getagiTAIR("/Users/gauravsablok/Desktop/release/test_ids.txt", \
                         "/Users/gauravsablok/Desktop/release/gene_association.tair")
{('AT5G27680', '501718175'),
 ('AT5G27680', '501741973'),
 ('AT5G27680', '501756966'),
 ('AT5G27680', '501780126')}                        
 help(getagiTAIR) # for detailed information
````

```python
from tairaccession.tairaccession import getagiDescription
getagiDescription(ids_file, association_file) 
#ids_file: path to the ids_file, 
#association_file: path to the association file
#provides the information about the association of the specific gene name and locus information on the agi
getagiDescription("/Users/gauravsablok/Desktop/release/test_ids.txt", \
                         "/Users/gauravsablok/Desktop/release/gene_association.tair")
{('AT5G27680', 'AT5G27680|AT5G27680.1|T1G16.10|T1G16_10'),
 ('AT5G27680', 'AT5G27680|AT5G27680.2'),
 ('AT5G27680', 'AT5G27680|AT5G27680.3'),
 ('AT5G27680', 'AT5G27680|AT5G27680.4'),
 ('AT5G27680', 'AT5G27680|AT5G27680.5'),
 ('AT5G27680', 'AT5G27680|AT5G27680.6'),
 ('AT5G27680', 'AT5G27680|RECQSIM|RECQ helicase SIM|T1G16.10|T1G16_10')}                         
help(getagiDescription) # for detailed information
````
```python
from tairaccession.functionalNames import functionalNames
# support for the phytozome gff files and functional name associations. also provides the option for the ids file to be specified.
functional = functionalNames
print(list(functional.keys())[:10])
['AT1G01010.1', 'AT1G01020.1', 'AT1G01020.2', 'AT1G01020.3', 'AT1G01020.4', 'AT1G01020.5', 'AT1G01020.6', 'AT1G01030.1', 'AT1G01030.2', 'AT1G01040.1']
print(list(functional.values())[:10])
['NAC domain containing protein 1', 'ARV1 family protein', 'ARV1 family protein', 'ARV1 family protein', 'ARV1 family protein', 'ARV1 family protein', 'ARV1 family protein', 'AP2/B3-like transcriptional factor family protein', 'AP2/B3-like transcriptional factor family protein', 'dicer-like 1']
# tairaccession has 2 inbuilt functions for the preparation of the update files as per the release of the phytozome.
```
```python
from tairaccession.geneNames import geneNames
genes = geneNames
print(list(genes.values())[:10])
['NAC001', 'ARV1', 'ARV1', 'ARV1', 'ARV1', 'ARV1', 'ARV1', 'NGA3', 'NGA3', 'DCL1']
print(list(genes.keys())[:10])
['AT1G01010.1', 'AT1G01020.1', 'AT1G01020.2', 'AT1G01020.3', 'AT1G01020.4', 'AT1G01020.5', 'AT1G01020.6', 'AT1G01030.1', 'AT1G01030.2', 'AT1G01040.1']
# a inbuilt gene name fetcher for any phytozome file. implemented # using the deque and provides faster iteration and prompt results. also provides the option for the ids file to be specified.
```
```python
from tairaccession.tairaccession import phytozomePacID
phytozomePacID("/Users/gauravsablok/Desktop/CodeTest/release/phytozome/Athaliana_447_Araport11.gene_exons.gff3", 
                   "/Users/gauravsablok/Desktop/CodeTest/release/all_ids.txt")
[('AT1G01010.1', 'pacid=37401853'),
 ('AT1G01020.1', 'pacid=37399351'),
 ('AT1G01020.4', 'pacid=37399352'),
 ('AT1G01020.3', 'pacid=37399353'),
 ('AT1G01020.5', 'pacid=37399354'),
 ('AT1G01020.2', 'pacid=37399355')]
help(phytozomePacID) # for detailed information
```
```python
import tairaccession
from tairaccession.geneNames import geneNames
GeneNames("/Users/gauravsablok/Desktop/CodeTest/release/all_ids.txt")
[('AT1G01010.1', 'NAC001'),
 ('AT1G01020.1', 'ARV1'),
 ('AT1G01020.2', 'ARV1'),
 ('AT1G01020.3', 'ARV1'),
 ('AT1G01020.4', 'ARV1'),
 ('AT1G01020.5', 'ARV1'),
 ('AT1G01020.6', 'ARV1'),
 ('AT1G01030.1', 'NGA3'),
 ('AT1G01030.2', 'NGA3')]
help(GeneNames) # for detailed documentation
```

```python
import tairaccession
from tairaccession.functionalNames import functionalNames
FunctionalNames("/Users/gauravsablok/Desktop/CodeTest/release/all_ids.txt")
[('AT1G01010.1', 'NAC domain containing protein 1'),
 ('AT1G01020.1', 'ARV1 family protein'),
 ('AT1G01020.2', 'ARV1 family protein'),
 ('AT1G01020.3', 'ARV1 family protein'),
 ('AT1G01020.4', 'ARV1 family protein'),
 ('AT1G01020.5', 'ARV1 family protein'),
 ('AT1G01020.6', 'ARV1 family protein')]
help(FunctionalNames) # for detailed documentation
```
```python
from tairaccession.tairaccession import visualizeAgiCDS
visualizeAgiCDS("/Users/gauravsablok/Desktop/CodeTest/release/gene_id.txt", 
        "/Users/gauravsablok/Desktop/CodeTest/release/TAIR10_GFF3_genes.gff", 
                                                    "arabidopsis", 10000, path)
help(visualizeAgiCDS) # for detail documentation
``````
![Visual](https://github.com/gauravcodepro/tairaccession/blob/main/gene_exon.png)

```python
from tairaccession.tairaccession import visualizeExons
visualizeAgiCDS("/Users/gauravsablok/Desktop/CodeTest/release/gene_id.txt", 
        "/Users/gauravsablok/Desktop/CodeTest/release/TAIR10_GFF3_genes.gff", 
                                                    "arabidopsis", 10000, path)
help(visualizeExons) # for detail documentation
``````
![Visual](https://github.com/gauravcodepro/tairaccession/blob/main/gene_exon.png)

```python
readTairNCBI("/Users/gauravsablok/Desktop/CodeCheck/release/ATH_GO_GOSLIM.txt", \
                  "/Users/gauravsablok/Desktop/CodeCheck/release/ATH_GO_GOSLIM.out", \
                          "2200935")
# this is just a reference id for checking the functionality
[('https://pubmed.ncbi.nlm.nih.gov/30356219/',
  'Nitrogen is an essential macronutrient for plant growth and basic metabolic processes. 
  The application of nitrogen-containing fertilizer increases yield, which has been a 
  substantial factor in the green revolution1. Ecologically, however, excessive application 
  of fertilizer has disastrous effects such as eutrophication2. A better understanding 
  of how plants regulate nitrogen metabolism is critical to increase plant yield and reduce 
  fertilizer overuse. Here we present a transcriptional regulatory network and twenty-one 
  transcription factors that regulate the architecture of root and shoot systems in response 
  to changes in nitrogen availability. Genetic perturbation of a subset of these transcription 
  factors revealed coordinate transcriptional regulation of enzymes involved in nitrogen metabolism. 
  Transcriptional regulators in the network are transcriptionally modified by feedback via 
  genetic perturbation of nitrogen metabolism. The network, genes and gene-regulatory modules 
  identified here will prove critical to increasing agricultural productivity.'),
 ('https://pubmed.ncbi.nlm.nih.gov/34562334/',
  'Unraveling gene function is pivotal to understanding the signaling cascades that control plant 
  development and stress responses. As experimental profiling is costly and labor intensive, 
  there is a clear need for high-confidence computational annotation. In contrast to detailed 
  gene-specific functional information, transcriptomics data are widely available for both model and 
  crop species. Here, we describe a novel automated function prediction method, which leverages 
  complementary information from multiple expression datasets by analyzing study-specific gene 
  co-expression networks. First, we benchmarked the prediction performance on recently characterized 
  Arabidopsis thaliana genes, and showed that our method outperforms state-of-the-art expression-based 
  approaches. Next, we predicted biological process annotations for known (n = 15 790) and unknown 
  (n = 11 865) genes in A. thaliana and validated our predictions using experimental protein-DNA and 
  protein-protein interaction data (covering >220 000 interactions in total), obtaining a set of 
  high-confidence functional annotations. Our method assigned at least one validated annotation to 
  5054 (42.6%) unknown genes, and at least one novel validated function to 3408 (53.0%) genes with 
  computational annotations only. These omics-supported functional annotations shed light on a 
  variety of developmental processes and molecular responses, such as flower and root development, 
  defense responses to fungi and bacteria, and phytohormone signaling, and help fill the information 
  gap on biological process annotations in Arabidopsis. An in-depth analysis of two context-specific 
  networks, modeling seed development and response to water deprivation, shows how previously 
  uncharacterized genes function within the respective networks. Moreover, our automated function 
  prediction approach can be applied in future studies to facilitate gene discovery for crop improvement.'),
 ('https://pubmed.ncbi.nlm.nih.gov/11118137/',
  'The completion of the Arabidopsis thaliana genome sequence allows a comparative analysis of 
  transcriptional regulators across the three eukaryotic kingdoms. Arabidopsis dedicates over 
  5% of its genome to code for more than 1500 transcription factors, about 45% of which are from 
  families specific to plants. Arabidopsis transcription factors that belong to families common 
  to all eukaryotes do not share significant similarity with those of the other kingdoms beyond the 
  conserved DNA binding domains, many of which have been arranged in combinations specific to each 
  lineage. The genome-wide comparison reveals the evolutionary generation of diversity in the 
  regulation of transcription.')]
```
## Contributing

Interested in contributing? Check out the contributing guidelines.  Please note that this project is released with a Code of Conduct.  By contributing to this project, you agree to abide by its terms.

## License

`tairaccession` was created by Gaurav Sablok, Universitat Potsdam, Germany. It is licensed under the terms of the MIT license.

Gaurav \
Academic Staff Member \
Bioinformatics \
Institute for Biochemistry and Biology \
University of Potsdam \
Potsdam,Germany

## Credits

`tairaccession` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) 
and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).
