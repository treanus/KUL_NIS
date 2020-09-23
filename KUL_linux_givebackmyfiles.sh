#!/bin/bash
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
#
# v0.1 - dd 23/09/2020

echo "This script will give back ownership of all fmriprep/mriqc directories. Type your password"
sudo chown -R $(id -u):$(id -g) mriqc* fmriprep* freesurfer*
echo "Done"

