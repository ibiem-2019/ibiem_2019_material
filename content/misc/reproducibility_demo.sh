#!/usr/bin/env bash

if ! [[ "$1" =~ ^(LOCAL|REMOTE|COMMAND)$ ]]; then
    printf "Run with LOCAL or REMOTE or COMMAND\n"
    printf "\tREMOTE: run RStudio\n"
    printf "\tLOCAL: run RStudio with access restricted to a web browser on the same machine\n"
    printf "\tCOMMAND: run demonstration dada2 and phyloseq pipeline, then exit\n"
    printf "\n\nMemory limitation can cause this to fail at the assignTaxonomy step.\n"
    printf "This can be a problem if there is not enough memory on the host machine,\n"
    printf "or if Docker is run with a memory cap. Docker Desktop on MacOS and Windows\n"
    printf "seems to default to a value that is too low"
    printf "\n\nAn optional second argument can be given to specify the base directory\n"
    exit 1
else
    REMOTE_OR_LOCAL=${1}
fi


BASE_DIR="/tmp/reproducible_demo_`date +%s`_tmp"
BASE_DIR=${2:-$BASE_DIR}

WORK_DIR="$BASE_DIR/work"
DATA_DIR="$BASE_DIR/data"
DOCKER_IMAGENAME="ibiem/docker_rstudio_ibiem2019:latest"
SEP_STRING="\n--------------------------------------------------\n"

# REMOTE_OR_LOCAL=$1
# MY_PASSWORD=${2:-"BadPa55word"}
MY_PASSWORD="`openssl rand -base64 16 | colrm 20`"

HOST_PORT=8787






mkdir -p $WORK_DIR $DATA_DIR

# git clone git@github.com:ibiem-2019/ibiem_2019_material.git ${WORK_DIR}/demo
git clone https://github.com/ibiem-2019/ibiem_2019_material.git ${WORK_DIR}/demo


printf "\n${SEP_STRING} Pulling docker image: $DOCKER_IMAGENAME ${SEP_STRING}"
# docker pull $DOCKER_IMAGENAME


if [ $REMOTE_OR_LOCAL == "COMMAND" ]; then
    printf "\n${SEP_STRING} STARTING Pipeline in Docker ${SEP_STRING}"
    docker run \
	   -it \
	   --rm \
	   -v ${WORK_DIR}:/home/guest \
	   -v ${DATA_DIR}:/data \
	   --user guest \
	   $DOCKER_IMAGENAME \
	   Rscript -e "rmarkdown::render('/home/guest/demo/content/lessons/run_everything.Rmd')"
    printf "\n${SEP_STRING} FINISHED Pipeline ${SEP_STRING}"
    printf "\n${SEP_STRING} Results output to ${WORK_DIR}/scratch"
    exit 0
elif [ $REMOTE_OR_LOCAL == "REMOTE" ]; then
    EXPOSE_PORT="${HOST_PORT}:8787"
    URL="http://`hostname -A | cut -f1 -d' '`:${HOST_PORT}/"
elif [ $REMOTE_OR_LOCAL == "LOCAL" ]; then
    EXPOSE_PORT="127.0.0.1:${HOST_PORT}:8787"
    URL="http://localhost:${HOST_PORT}/"
else
    printf "Run with LOCAL or REMOTE or COMMAND\n"
    printf "\tREMOTE: run RStudio\n"
    printf "\tLOCAL: run RStudio with access restricted to a web browser on the same machine\n"
    printf "\tCOMMAND: run demonstration dada2 and phyloseq pipeline, then exit\n"
    exit 1
fi
echo "expose argument: ${EXPOSE_PORT}"


printf "\n${SEP_STRING} Running Docker ${SEP_STRING}"
DOCKER_CMD="docker run --name ibiem -d -p ${EXPOSE_PORT} -e USERPASS=${MY_PASSWORD} -v ${WORK_DIR}:/home/guest -v ${DATA_DIR}:/data $DOCKER_IMAGENAME"
echo $DOCKER_CMD
$DOCKER_CMD

printf "${SEP_STRING}URL:\t${URL}\nUsername: guest\nPassword: ${MY_PASSWORD} ${SEP_STRING}"
printf "${SEP_STRING}In RStudio:\n1. New Project: https://github.com/ibiem-2019/ibiem_2019_material.git"


printf "\n${SEP_STRING} To clean up: \n\n"
echo "rm -rf $BASE_DIR"
echo "docker rm -f \`docker ps -aq\`"
echo "docker rmi \`docker images -aq\`"
printf "${SEP_STRING}\n"

