#Installation

We recommed you install the limix based QTL mapping pipeline in a seperate conda enviroment.

To do so please start by making a conda enviroment, [conda](https://conda.io/docs/index.html):.

`conda create -n limix_qtl python=2.7`
`source activate limix_qtl`
`conda install -c conda-forge limix`
`pip install tables`
`conda install -c conda-forge pandas-plink`
`pip install limix --upgrade`


//Doesn't work yet.
`git clone https://github.com/PMBio/hipsci_pipeline.git`


To create a limix 1.1 enviroment

`conda create -n limix_1_1_qtl python=2.7`
`source activate limix_1_1_qtl`
`git clone https://github.com/limix/limix`
`cd limix`
`git checkout release/1.1.0`
`conda install -c conda-forge liknorm`
`pip install liknorm-py pandas-plink glimix-core limix-core limix-legacy optimix --upgrade --upgrade-strategy='only-if-needed'`
`python setup.py develop`