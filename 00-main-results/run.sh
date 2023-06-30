#!/usr/bin/env bash
set -ex

# If not running in a CI, get commit ID from git
if [ -z "${CI_COMMIT_SHORT_SHA}" ]
then
    COMMIT=$(git rev-parse --short HEAD)
else
    COMMIT=${CI_COMMIT_SHORT_SHA}
fi

BUILD_DATE=$(date -u +'%Y-%m-%dT%H:%M:%SZ')
BUILD_DATE_INT=$(date -u +'%Y%m%d')
IMAGE_NAME="snapanalysis"

CONTAINER_NAME="snapanalysis_cont_${BUILD_DATE_INT}_${COMMIT}"

# https://stackoverflow.com/questions/36335186/bash-exit-and-cleanup-on-error
clean_up () {
    ARG=$?
    if [ $ARG -ne 0 ]
    then
        echo "Got error, code=$ARG, cleaning up"
        docker stop ${CONTAINER_NAME} || echo "Tried stopping ${CONTAINER_NAME}, didn't work"
        docker rm ${CONTAINER_NAME} || echo "Tried removing ${CONTAINER_NAME}, didn't work"
        exit $ARG
    fi
}
trap clean_up EXIT

docker build --build-arg=COMMIT=${COMMIT}  --build-arg=BUILD_DATE=${BUILD_DATE} -f Dockerfile -t ${IMAGE_NAME}:${BUILD_DATE_INT} -t ${IMAGE_NAME}:${COMMIT} -t ${IMAGE_NAME}:latest .
docker create -it --name ${CONTAINER_NAME} snapanalysis:${COMMIT}
docker cp ./sort_output.sh ${CONTAINER_NAME}:/root/sort_output.sh
docker start ${CONTAINER_NAME}
docker exec ${CONTAINER_NAME} /root/sort_output.sh /output/ /output/sorted/ y
docker stop ${CONTAINER_NAME}
docker cp ${CONTAINER_NAME}:/output/. output
docker cp ${CONTAINER_NAME}:/data/. output/data-cache/
docker rm ${CONTAINER_NAME}

# Compress data for web into archive
echo "$BUILD_DATE\n$COMMIT" > output/web/VERSION
tar -czvf output/marcs.for-web.tar.gz output/web

echo "$BUILD_DATE\n$COMMIT" > output/last_successful_build.txt
