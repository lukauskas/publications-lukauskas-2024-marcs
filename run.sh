#!/usr/bin/env bash
set -ex

# https://stackoverflow.com/questions/36335186/bash-exit-and-cleanup-on-error
clean_up () {
    ARG=$?
    echo "Got error, cleaning up"
    docker stop snapanalysis_cont || echo "Tried stopping snapanalysis_cont, didn't work"
    docker rm snapanalysis_cont || echo "Tried removing snapanalysis_cont, didn't work"
    exit $ARG
}
trap clean_up EXIT

docker build --build-arg=COMMIT=$(git rev-parse --short HEAD)  --build-arg=BUILD_DATE=$(date -u +'%Y-%m-%dT%H:%M:%SZ') -f Dockerfile -t snapanalysis .
docker create -it --name snapanalysis_cont snapanalysis
docker cp ./sort_output.sh snapanalysis_cont:/root/sort_output.sh
docker start snapanalysis_cont
docker exec snapanalysis_cont /root/sort_output.sh /output/ /output/sorted/ y
docker stop snapanalysis_cont
docker cp snapanalysis_cont:/output/. output
docker rm snapanalysis_cont
