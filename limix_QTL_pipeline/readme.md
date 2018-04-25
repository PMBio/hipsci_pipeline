#Installation

We recommed you install the limix based QTL mapping pipeline in a seperate conda enviroment.

//Doesn't work yet.
`git clone https://github.com/PMBio/hipsci_pipeline.git`


To create a limix development enviroment

`conda create -n limix_qtl_develop python=3 anaconda`

`source activate limix_qtl_develop`

`bash <(curl -fsSL https://raw.githubusercontent.com/limix/limix/develop/install)`

`conda install -c anaconda pytest`

`conda install pytables`


NB. be sure to be in a folder where you can download files to and there is no folder called limix.

Some filesystems have locking disabled to write to be able to use the tool use:
export HDF5_USE_FILE_LOCKING=FALSE
