#!/bin/bash
# @ Stefan Sunaert & Ahmed Radwan- UZ/KUL - stefan.sunaert@uzleuven.be
what_to_build=$1

if [ "$what_to_build" = "" ]; then

    echo "Use KUL_build_singularity what_to_build "
    echo "  what to build could be e.g. fmriprep:latest or mriqc:0.12.4"
    exit 0

fi

cwd=$(pwd)
docker run --privileged -t --rm \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v ${cwd}:/output \
    singularityware/docker2singularity \
    poldracklab/fmriprep:latest
