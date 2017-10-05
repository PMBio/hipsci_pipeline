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

`git clone https://github.com/limix/limix`

`cd limix`

`git checkout develop`

`conda install -c conda-forge liknorm bgen asciitree hcephes humanfriendly tabulate`

`conda install libffi statsmodels`

`conda install -c anaconda joblib`

`pip install liknorm-py pandas-plink glimix-core limix-core limix-legacy optimix tables --upgrade --upgrade-strategy='only-if-needed'`

`python setup.py develop`
