#MicrobiomeGS2 - constrained based community metabolic modeling in R

MicrobiomeGS2 uses the approach of steadycom and ports it to methods implemented in R. 

## Installation

### Using Conda


Install the conda dependencies
`curl https://raw.githubusercontent.com/Waschina/MicrobiomeGS2/master/requirements.yml
conda env create -f requirements.yml
rm requirements.yml
conda activate MicrobiomeGS2`

Install missing R packages
install.packages(c("https://cran.r-project.org/src/contrib/Archive/sybil/sybil_2.2.0.tar.gz", "https:/,"/cran.r-project.org/src/contrib/Archive/sybilSBML/sybilSBML_3.1.2.tar.gz"), repos = NULL)

`Rscript -e "install.packages(c( 'https://cran.r-project.org/src/contrib/Archive/sybil/sybil_2.2.0.tar.gz', 'https://cran.r-project.org/src/contrib/Archive/sybilSBML/sybilSBML_3.1.2.tar.gz'), repos = NULL)"
Rscript -e "`devtools::install_github('Waschina/MicrobiomeGS2')"`

