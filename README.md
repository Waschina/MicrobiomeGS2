# MicrobiomeGS2 

***constraint-based community metabolic modeling in R***

MicrobiomeGS2 uses the approach of steadycom and ports it to methods implemented in R. 

## Installation

### Using Conda


Install the conda dependencies within a new environment with the name `MicrobiomeGS2`.

```sh
wget https://raw.githubusercontent.com/Waschina/MicrobiomeGS2/master/requirements.yml
conda env create -f requirements.yml
rm requirements.yml
conda activate MicrobiomeGS2
```

Install the MicrobiomeGS2 package and dependencies required for SBML support

```sh
Rscript -e "devtools::install_github('Waschina/MicrobiomeGS2')"
wget https://cran.r-project.org/src/contrib/Archive/sybilSBML/sybilSBML_3.1.2.tar.gz
R CMD INSTALL --configure-args=" \
--with-sbml-include=$CONDA_PREFIX/include \
--with-sbml-lib=$CONDA_PREFIX/lib" sybilSBML_3.1.2.tar.gz
rm sybilSBML_3.1.2.tar.gz
```

## CPLEX support

Download and install cplex from the [IBM homepage](https://www.ibm.com/support/pages/downloading-ibm-ilog-cplex-optimization-studio-v12100). **Academic users** can get a free licence for CPLEX and [download the software for free](https://github.com/academic-initiative/documentation/blob/main/academic-initiative/how-to/How-to-download-IBM-ILOG-CPLEX/readme.md). Get the cplexAPI package and link the library.

```sh
wget https://cran.r-project.org/src/contrib/Archive/cplexAPI/cplexAPI_1.4.0.tar.gz
R CMD INSTALL --configure-args=" \
   --with-cplex-dir=/path/to/cplex12.10/CPLEX_Studio/cplex" \
    cplexAPI_1.4.0.tar.gz
rm cplexAPI_1.4.0.tar.gz
```
