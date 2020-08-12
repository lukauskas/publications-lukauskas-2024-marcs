#!/usr/bin/env bash
set -e
docker build --build-arg=COMMIT=$(git rev-parse --short HEAD)  --build-arg=BUILD_DATE=$(date -u +'%Y-%m-%dT%H:%M:%SZ') -f Dockerfile -t snapanalysis .
docker create -it --name snapanalysis_cont snapanalysis
docker cp snapanalysis_cont:/output/. output
docker rm snapanalysis_cont
