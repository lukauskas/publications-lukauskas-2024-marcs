MARCS Analysis Code
==============================

This repository hosts the code accompanying the
Modification Atlas of Regulation by Chromatin States (MARCS)

The main analysis steps and figures can be generated from this pipeline.

## Prerequisites

### Operating system

The analysis can only be built on Linux and Mac OS X.

### Docker

Docker daemon has to be installed on the system before starting.

* MAC OS X:  [follow these instructions](https://docs.docker.com/docker-for-mac/install/)
    * Note that Docker for Mac OS X comes with very restrictive default settings for resource usage.
      You need to increase the available memory to at least 8GB (Preferences > Resources > Memory 8GB) to run this pipeline.
      Additionally, increase the number of available CPU cores for faster runtimes.
* linux: use your distro's package manager.


### Helvetica

In order to build this analysis, especially the figures,
you need to source a copy of Helvetica font.
As it is proprietary it is not included in this repository.

The simplest way to obtain it is from your Mac OS X.

An abridged version of [Olga Botvinnik's](https://olgabotvinnik.com/blog/how-to-set-helvetica-as-the-default-sans-serif-font-in/) and [gree2](https://gree2.github.io/python/2015/04/27/python-change-matplotlib-font-on-mac) guides is below:

1. Find `Helvetica` on your system. It will either be:
    * `/System/Library/Fonts/Helvetica.ttc` (Mojave+)
    * or `/System/Library/Fonts/Helvetica.dfont` (older?)
2. Download [`dfontsplitter`](http://peter.upfold.org.uk/projects/dfontsplitter) (`ttc` or `dfont` files) or
  `fondu` (`dfont` files)
3. Use `dfontsplitter` or `fondu` to convert `Helvetica` into `ttf` format.
4. Put the `ttf` files in the `fonts` directory

The Dockerfile will do the rest for you.

## Structure

The directory tree is split into three branches

1. `data` - stores raw data as well as extra downloaded resources
2. `src` - stores the main source code for analyses
3. `extra-figures-notebooks` - stores jupyter notebooks for most figures used in the paper

### `data` directory

The starting point of this analysis is the `raw` subdirectory of the data in which the output of `MaxQuant` is stored in the `proteinGroups.txt.gz` file. Additionally, the `raw` directory hosts metadata about the experiments as well as manually curated list of protein complexes that are present in the data.

In the `raw/manual-overrides` subdirectory, manual overrides of network coordinates are stored for visualisation. As the network layout procedures are stochastic, the manual overrides are there to ensure the produced network looks the same every time.

In `downloaded` subdirectory the data from orthogonal sources is stored. Namely, this directory hosts data from [EpiFactors](https://epifactors.autosome.ru/), [EBI Complex Portal](https://www.ebi.ac.uk/complexportal/), [BioGRID](https://thebiogrid.org/) and [MyGene.info](https://mygene.info/).
In theory, the `Makefile` in this directory will re-download this data. In practice, the data location and format will most likely be out of date and the scripts will therefore need to be adjusted.

### `src` directory

The source directory hosts the main python package for analysis
which is called `snapanalysis` (Silac Nucleosome Affinity Pulldown (SNAP) Analysis). The package can be installed with `setup.py` script, but be wary about requirements. The requirements in docker container are most up-to-date.

The module is divided into parts as follows:

* `snapanalysis.preprocessing` contains the scripts that parse data in `raw` folder sit in the preprocessing directory. This is where unique human readable labels (`Gene label` column) get assigned to proteins.
* `snapanalysis.external` contains the scripts that parse external datasets.
* `snapanalysis.models` holds the main analyses in the paper with:
    * `snapanalysis.models.enrichment` is the enrichment decomposition step which corrects for H/L mixing.
    * `snapanalysis.models.network` is the network analysis submodule
    * `snapanalysis.models.ptm_response` is the LIMMA model for response to chromatin modifications.
* `snapanalysis.visualisation` holds the visualisation scripts that are frequently used such as the heatmap visualisations with nucleosome composition header.
* `snapanalysis.helpers` was intended to be a general helper module but ended up only having a single function in it

While this model can be installed without the docker, it is not advisable to run the analyses without the corresponding docker environment as it is very tightly integrated into it.

### `extra-figures-notebooks`

This directory hosts the `ipynb` notebooks for figures used in the paper. The directories are split into topics.

The notebooks assume the pipeline has already been run and the results are available.

## Running the pipeline

The pipeline is encoded into the `Dockerfile`.
The scripts in the dockerfile will run the analyses in correct order, will load `extra-figures-notebooks` and run them as well.

The output is then stored in `/output` in the container and can be downloaded this way.

A helper script `run.sh` automates the build steps and saves output to `output` directory.

To run the whole pipeline, therefore, it is sufficient to just type:

```
./run.sh
```

Default name for the built image is `snapanalysis`.
This image will be saved locally and is pre-requisite to build the web interface.

## Development

There are two ways to run the pipeline in development mode, each serving a different purpose.
Quick proof-of-concept analyses and creation of figures are best done in `jupyter` environment.

Stable models, on the other hand, and other major resources should be incorporated into `snapanalysis` python package.
This is best done in PyCharm environment.


### Development in Jupyter

Jupyter notebooks are stored in `extra-figures-notebooks` directory.
To hot-load them use the script `run_dev_notebooks.sh`:

```
`./run_dev_notebooks.sh
```

Notebooks will be accessible on [http://localhost:8889](http://localhost:8889).
Make sure to copy the token that is printed into your terminal.

Note that this script will enter the last *successfully* built container.
If you haven't run `./run.sh` for a while, or the last few runs crashed,
the notebooks will become out-of-sync.

Since the notebooks are loaded as volumes, changes you make will persist on disk.

### Development in PyCharm

PyCharm can be set-up to work with the main `snapanalysis` package.
It is best to create a python3 virtual environment exclusively for this project. The package should then be installed into this environment in editable mode by navigating to `src` directory and typing `pip install -e .`. There might be a few conflicts regarding the dependancies which will have to be solved then.

Inside the PyCharm, setting `src` directory as "library root" will then set the paths up correctly and the code autocompletion and search should start to work.

Note that jupyter notebooks cannot be developed from PyCharm, but the main library (`snapanalysis` package) can, and vica-versa.

It should be possible to setup PyCharm to run inside the Docker environment, but I never found the time to do that, but the scripts provide a workaround.

## Other links

The mass spectrometry proteomics data used in this analysis
have been deposited to the ProteomeXchange Consortium via the PRIDE (Perez-Riverol, 2019) partner
repository with the dataset identifier [PXD018966](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD018966).

## References

Perez-Riverol Y, Csordas A, Bai J, Bernal-Llinares M, Hewapathirana S, Kundu DJ, Inuganti A, Griss J, Mayer G, Eisenacher M, Pérez E, Uszkoreit J, Pfeuffer J, Sachsenberg T, Yilmaz S, Tiwary S, Cox J, Audain E, Walzer M, Jarnuczak AF, Ternent T, Brazma A, Vizcaíno JA (2019). The PRIDE database and related tools and resources in 2019: improving support for quantification data. Nucleic Acids Res 47(D1):D442-D450 (PubMed ID: 30395289).
