#Installation

We recommed you install the limix based QTL mapping pipeline in a seperate conda enviroment.

To do so please start by making a conda enviroment, [conda](https://conda.io/docs/index.html):.

`conda create -n limix_qtl python=2.7`

`source activate limix_qtl`

`conda install -c conda-forge limix bgen pandas-plink`

`pip install tables`

`pip install limix --upgrade`



//Doesn't work yet.
`git clone https://github.com/PMBio/hipsci_pipeline.git`


To create a limix development enviroment

`conda create -n limix_qtl_develop python=3 anaconda`

`source activate limix_qtl_develop`

`bash <(curl -fsSL https://raw.githubusercontent.com/limix/limix/develop/install)`

`conda install -c anaconda pytest`

`conda install -c ska tables`


NB. be sure to be in a folder where you can download files to and there is no folder called limix.
