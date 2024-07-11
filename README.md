# tairaccession

<div align = "justify">
  
- tairaccession: This allows you to analyze the tair accession ids easily and you can see the sample datasets for each of the import function with in the tests directory and the corresponding format are available from [TAIR](https://www.arabidopsis.org). This package has added utilities for analyzing also the Phytozome datasets and also you can plot the desried genes of interest.
- You can also analyze [Phytozome Araport](https://phytozome-next.jgi.doe.gov) also using this package.
- The package is under release at PyPI package repository.
- I have updated this package with additional support for the phytozome in addition to the tair.
- In an update to this package, few functions on plotting the coding regions, gene regions and exons have been added.
- There are additional functions such ```prepareFunctionalNamePhytozome```and ```preparegeneNamePhytozome```which will automatically prepare the files as per the release of the phytozome and tair.
- The web documentation is located at [tairaccession](https://gauravcodepro.github.io/tairaccession)

If you have any questions, please contact at gauravcodepro@gmail.com.
## Installation
```bash
$ pip install tairaccession
import tairaccession
print(tairaccession.__version__)
```

## License
`tairaccession` was created by Gaurav Sablok, Universitat Potsdam, Germany. It is licensed under the terms of the MIT license.

Gaurav 

## Credits
`tairaccession` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) 
and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).
