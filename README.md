#MicrobiomeGS2 - constrained based community metabolic modeling in R

MicrobiomeGS2 uses the approach of steadycom and ports it to methods implemented in R. 

## Installation

### Using Conda


Install the conda dependencies

```
curl https://raw.githubusercontent.com/Waschina/MicrobiomeGS2/master/requirements.yml
conda env create -f requirements.yml
rm requirements.yml
conda activate MicrobiomeGS2
```

Install missing R packages

```
Rscript -e "install.packages(c( 'https://cran.r-project.org/src/contrib/Archive/sybil/sybil_2.2.0.tar.gz', 'https://cran.r-project.org/src/contrib/Archive/sybilSBML/sybilSBML_3.1.2.tar.gz'), repos = NULL)"
Rscript -e "devtools::install_github('Waschina/MicrobiomeGS2')"
```

## CPLEX support

Download and install cplex from the [IBM homepage](https://www.ibm.com/support/pages/downloading-ibm-ilog-cplex-optimization-studio-v12100). Get the cplexAPI package and link the library.

```
wget https://cran.r-project.org/src/contrib/Archive/cplexAPI/cplexAPI_1.4.0.tar.gz
R CMD INSTALL --configure-args=" \
   --with-cplex-dir=/path/to/cplex12.10/CPLEX_Studio/cplex" \
    cplexAPI_1.4.0.tar.gz
rm cplexAPI_1.4.0.tar.gz
```
