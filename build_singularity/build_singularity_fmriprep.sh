
cwd=$(pwd)
docker run --privileged -t --rm \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v /absolute/path/to/output/folder:/output \
    singularityware/docker2singularity \
    poldracklab/fmriprep:<version>
