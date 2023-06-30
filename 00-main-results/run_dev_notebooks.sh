#!/usr/bin/env bash
set -e
docker run --rm -w /notebooks -it -p 8889:8889 --name snapanalysis_notebooks_cont -v $(pwd)/src:/snap:ro -v $(pwd)/extra-figures-notebooks:/notebooks -v $(pwd)/output:/output snapanalysis jupyter notebook --port 8889 --ip=0.0.0.0 --allow-root
