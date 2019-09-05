#!/usr/bin/env bash

DEMO_BASE="/tmp/reproducible_demo"
WORK_DIR="$DEMO_BASE/work"
DATA_DIR="$DEMO_BASE/data"
DOCKER_IMAGENAME="ibiem/docker_rstudio_ibiem2019"


mkdir -p $WORK_DIR $DATA_DIR

SEP_STRING="\n--------------------------------------------------\n"

printf "\n${SEP_STRING} Pulling docker image: $DOCKER_IMAGENAME ${SEP_STRING}"
docker pull $DOCKER_IMAGENAME

printf "\n${SEP_STRING} Running Docker ${SEP_STRING}"
docker run \
    -i \
    -v ${WORK_DIR}:/home/guest \
    -v ${DATA_DIR}:/data \
    $DOCKER_IMAGENAME head /home/guest/demo/content/lessons/run_everything.Rmd


# Rscript -e "rmarkdown::render('output_file.Rmd')"
printf "\n${SEP_STRING} Running Docker ${SEP_STRING}"
docker run \
    -i \
    -v ${WORK_DIR}:/home/guest \
    -v ${DATA_DIR}:/data \
    $DOCKER_IMAGENAME Rscript -e "rmarkdown::render('/home/guest/demo/content/lessons/run_everything.Rmd')"

