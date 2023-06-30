# MARCS integration with H3K4me1/3 ChIP-MS

See [report.md](report.md) for information, as well as the `ipynb` notebooks (in the current directory), the [outputs](outputs) and the [data](data) directory (inputs and docs).

## Requirements

The code can be run on any OSX/unix machine with a functioning conda installation.

Set up the conda environment:

```
conda env create -f environment.yaml
```

Activate the conda environment:

```
conda activate marcs-chip-ms-h3k4me1-h3k4me3
```

Run the pipeline:

```
./run.sh
```

(or edit notebooks manually using `jupyter lab`)

For software versions see inside the environment.yaml file
