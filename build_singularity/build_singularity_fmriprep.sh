#!/bin/bash
# @ Stefan Sunaert & Ahmed Radwan- UZ/KUL - stefan.sunaert@uzleuven.be

cwd=$(pwd)
docker run --privileged -t --rm \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v ${cwd}:/output \
    singularityware/docker2singularity \
    poldracklab/fmriprep:latest
