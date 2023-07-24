# tairaccession

tairaccession: This package has been developed at University of Pretoria, South Africa and licensed to University of Pretoria, South Africa. This allows you to analyze the tair accession ids easily and you can see the sample datasets for each of the import function with in the tests directory and the corresponding format are available from [TAIR](https://www.arabidopsis.org). This package has added utilities for analyzing also the Phytozome datasets and also you can plot the desried genes of interest. You can also analyze Phytozome Araport also using this package.
If you have any questions, please contact at sablokg@gmail.com. The package is under release at PyPi package repository. 

## Installation

```bash
$ pip install tairaccession
import tairaccession
print(tairaccession.__version__)
```

## Usage

`tairaccession` can be used to access the tairIDs and also used for the conversion of the tair data and obtaining the relative information
from the tair accessions aka AGIs specificed in a file. The tair accession has the following options:


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
 ('A0A1P8BBY0', ['locus:2180255', 'AT5G27680']),
 ('A0A1P8BBY3', ['locus:2180255', 'AT5G27680']),
 ('A0A654F550', ['locus:2096444', 'AT3G03550']),
 ('Q0WS90', ['locus:2096444', 'AT3G03550']),
 ('A0A654EWS4', ['locus:2010796', 'AT1G64990']),
 ('A0A654F9L3', ['locus:2091623', 'AT3G22240'])]  
help(uniprotTair) # for detailed information
````
```python
from tairaccession.tairaccession import agiCoordinates
agiCoordinates(ids_file, gff_file, gene_type = None)  
#ids_file: path to the ids_file #gff_file: path to the gff file 
#gene_type: this is a #keyworded argument to the agi_coordinates. 
#if gene_type is gene: it will return the gene coordinates for the agi in the file, 
#if gene_type #is exon: it will return the exon coordinates for the agi in the file,
#if gene_type is three_prime_UTR: it will return the three_prime_UTR 
#if gene_type is five_prime_UTR: it will return the five_prime_UTR coordinates for the agi in the file,
agiCoordinates("/Users/gauravsablok/Desktop/release/test_ids.txt", \
                    "/Users/gauravsablok/Desktop/release/TAIR10_GFF3_genes.gff", \
                                                                gene_type = "gene")
[[['AT1G64990'], [24139907, 24146189]],
 [['AT3G03550'], [850108, 851590]],
 [['AT3G22240'], [7863616, 7864819]],
 [['AT3G24715'], [9025849, 9029948]],
 [['AT5G27680'], [9793790, 9798961]]]

agiCoordinates("/Users/gauravsablok/Desktop/release/test_ids.txt", \
                    "/Users/gauravsablok/Desktop/release/TAIR10_GFF3_genes.gff", \
                                                                  gene_type = "exon")
[[['AT1G64990', 24145652, 24145980]],
 [['AT1G64990', 24145055, 24145155]],
 [['AT1G64990', 24144562, 24144669]],
 [['AT1G64990', 24143921, 24144041]],
 [['AT1G64990', 24143323, 24143452]],
 [['AT1G64990', 24142736, 24142814]],
 [['AT1G64990', 24142038, 24142126]],
 [['AT1G64990', 24141775, 24141867]],
 [['AT1G64990', 24141637, 24141663]],
 [['AT1G64990', 24141437, 24141523]],
 [['AT1G64990', 24141031, 24141102]],
 [['AT1G64990', 24140542, 24140672]],
 [['AT1G64990', 24139989, 24140336]],
 [['AT1G64990', 24146019, 24146189]],
 [['AT1G64990', 24145652, 24145915]],
 [['AT1G64990', 24145055, 24145155]],
 [['AT1G64990', 24144562, 24144669]],
 [['AT1G64990', 24143921, 24144041]],
 [['AT1G64990', 24143323, 24143452]],
 [['AT1G64990', 24142736, 24142814]],
 [['AT1G64990', 24142038, 24142126]],
 [['AT1G64990', 24141775, 24141867]],
 [['AT1G64990', 24141637, 24141663]],
 [['AT1G64990', 24141437, 24141523]],
 [['AT1G64990', 24141031, 24141102]],
 [['AT1G64990', 24140542, 24140672]],
 [['AT1G64990', 24139907, 24140336]]]
 agiCoordinates("/Users/gauravsablok/Desktop/release/test_ids.txt", \
                     "/Users/gauravsablok/Desktop/release/TAIR10_GFF3_genes.gff", \
                                                       gene_type = "three_prime_UTR")
[['AT1G64990', 24139989, 24140182],
 ['AT1G64990', 24139907, 24140182]] 
agiCoordinates("/Users/gauravsablok/Desktop/release/test_ids.txt", \
                     "/Users/gauravsablok/Desktop/release/TAIR10_GFF3_genes.gff", \
                                                       gene_type = "five_prime_UTR") 
[[['AT1G64990', 24145867, 24145980]],
 [['AT1G64990', 24146019, 24146189]],
 [['AT1G64990', 24145867, 24145915]]]
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
## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`tairaccession` was created by Gaurav Sablok, University of Pretoria, South Africa. It is licensed under the terms of the MIT license.

## Credits

`tairaccession` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).
