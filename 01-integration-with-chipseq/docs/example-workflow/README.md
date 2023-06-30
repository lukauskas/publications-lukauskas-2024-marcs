# An example of Mutual Information Statistic Calculation

This directory serves as a simplified example of the processing that is done by the pipeline, so it is a bit easier to understand.
The outputs of this script, for two BED files should be identical to those of the pipeline, and this is checked in the last step.

See the [ExampleProcessing notebook](ExampleProcessing.ipynb) for more info.

## To run this example

Install and activate the conda environment as usual:

```
conda env create -f environment.yaml
```

```
conda env activate marcs-chipseq-example
```

Then you should be able to access the jupyter notebook through 

```
jupyter lab
```

The last step of the example notebook assumes that the pipeline has been run and the outputs were generated.
If you want to run these cells, make sure the outputs are generated and can be accessed by the script.

