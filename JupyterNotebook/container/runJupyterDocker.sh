#!/bin/bash
# A script to run the docker image with the jupyter notebook

# The path to the local directory containing the data and the pipeline
DATA_DIR=_____
PIPELINE_DIR= _____/pipeline
UTILS_DIR=_____/utils
# The path to the working directory in the container.
# ** WARNING: 
#    It is best to leave this path unchanged, so that it matches the path for the
#    working directory in the Dockerfile that was used to build the container.
WORK_DIR=/tmp/work

# The name and version of the container to run
CONTAINER_NAME=analysis-pipeline
# If a container version is not specified, version 1.1.0 will be used
if [ -z "$1" ]
then
    CONTAINER_VERSION=1.1.0
else
    CONTAINER_VERSION=$1
fi

docker run -it --rm -u root -e GRANT_SUDO=yes -p 8888:8888 -p 4545:4545 -v $PIPELINE_DIR:$WORK_DIR/pipeline -v $DATA_DIR:$WORK_DIR/data  -v $UTILS_DIR:$WORK_DIR/utils $CONTAINER_NAME:$CONTAINER_VERSION
