#!/usr/bin/env bash
set -e

OUTPUT_DIR='./outputs'

if [ ! -d "${OUTPUT_DIR}" ]; then
  mkdir -p ${OUTPUT_DIR}
else
    echo "Note that output directory ${OUTPUT_DIR} already exists, we will overwrite, but not delete outputs"
fi



jupyter nbconvert --output-dir="${OUTPUT_DIR}/notebooks" --ExecutePreprocessor.timeout=600 --to notebook --execute "01-extracting.ipynb"
jupyter nbconvert --output-dir="${OUTPUT_DIR}/notebooks" --ExecutePreprocessor.timeout=600 --to notebook --execute "02-linking-to-MARCS.ipynb"
jupyter nbconvert --output-dir="${OUTPUT_DIR}/notebooks" --ExecutePreprocessor.timeout=600 --to notebook --execute "03-transformation-and-modelling.ipynb"
jupyter nbconvert --output-dir="${OUTPUT_DIR}/notebooks" --ExecutePreprocessor.timeout=600 --to notebook --execute "04-integrating-results-with-MARCS.ipynb"
jupyter nbconvert --output-dir="${OUTPUT_DIR}/notebooks" --ExecutePreprocessor.timeout=600 --to notebook --execute "05-excel-output.ipynb"
