#!/usr/bin/env bash
set -e

OUTPUT_DIR='./outputs'

if [ ! -d "${OUTPUT_DIR}" ]; then
  mkdir -p ${OUTPUT_DIR}
else
    echo "Note that output directory ${OUTPUT_DIR} already exists, we will overwrite, but not delete outputs"
fi

# For some reason this is needed
ulimit -n 4096

jupyter nbconvert --output-dir="${OUTPUT_DIR}/notebooks" --ExecutePreprocessor.timeout=600 --to notebook --execute "01-extracting.ipynb"
papermill "02-transformation-and-modelling.ipynb" "${OUTPUT_DIR}/notebooks/02-transformation-and-modelling-long-linkers-prom.ipynb" -p DATASET long-linkers-prom
papermill "02-transformation-and-modelling.ipynb" "${OUTPUT_DIR}/notebooks/02-transformation-and-modelling-long-linkers-enh.ipynb" -p DATASET long-linkers-enh
papermill "02-transformation-and-modelling.ipynb" "${OUTPUT_DIR}/notebooks/02-transformation-and-modelling-short-linkers.ipynb" -p DATASET short-linkers

papermill "03-excel-output.ipynb" "${OUTPUT_DIR}/notebooks/03-excel-output-long-linkers-enh.ipynb" -p DATASET long-linkers-enh
papermill "03-excel-output.ipynb" "${OUTPUT_DIR}/notebooks/03-excel-output-long-linkers-prom.ipynb" -p DATASET long-linkers-prom
papermill "03-excel-output.ipynb" "${OUTPUT_DIR}/notebooks/03-excel-output-short-linkers.ipynb" -p DATASET short-linkers

papermill "04-plots.ipynb" "${OUTPUT_DIR}/notebooks/04-plots-long-linkers-enh.ipynb" -p DATASET long-linkers-enh
papermill "04-plots.ipynb" "${OUTPUT_DIR}/notebooks/04-plots-long-linkers-prom.ipynb" -p DATASET long-linkers-prom
papermill "04-plots.ipynb" "${OUTPUT_DIR}/notebooks/04-plots-short-linkers.ipynb" -p DATASET short-linkers
