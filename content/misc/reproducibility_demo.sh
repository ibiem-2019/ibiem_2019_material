DEMO_BASE="/tmp/reproducible_demo"
WORK_DIR="$DEMO_BASE/work"
DATA_DIR="$DEMO_BASE/data"
DOCKER_IMAGENAME="ibiem/docker_rstudio_ibiem2019"
SEP_STRING="\n--------------------------------------------------\n"


REMOTE_OR_LOCAL=${1:-"LOCAL"}
MY_PASSWORD=${2:-"BadPa55word"}

HOST_PORT=8787

if [ $REMOTE_OR_LOCAL == "REMOTE" ]; then
    EXPOSE_PORT="${HOST_PORT}:8787"
    URL="http://`hostname -A | cut -f1 -d' '`:${HOST_PORT}/"
else
    EXPOSE_PORT="127.0.0.1\:${HOST_PORT}\:8787"
    URL="http://localhost:${HOST_PORT}/"
fi
echo "expose argument: ${EXPOSE_PORT}"




mkdir -p $WORK_DIR $DATA_DIR

git clone git@github.com:ibiem-2019/ibiem_2019_material.git ${WORK_DIR}/demo
# git clone https://github.com/ibiem-2019/ibiem_2019_material.git ${WORK_DIR}/demo


printf "\n${SEP_STRING} Pulling docker image: $DOCKER_IMAGENAME ${SEP_STRING}"
docker pull $DOCKER_IMAGENAME

# 
# printf "\n${SEP_STRING} Downloading data from DDS ${SEP_STRING}"
# 
# ddsclient download -p HTS_course --include hts_2019_data/hts2019_pilot_rawdata/21_2019_P_M1_S21_L002_R1_001.fastq.gz $DATA_DIR
# 
# 
# printf "\n${SEP_STRING} Cloning repo from gitlab ${SEP_STRING}"
# 


printf "\n${SEP_STRING} Running Docker ${SEP_STRING}"
DOCKER_CMD="docker run --name ibiem -d -p ${EXPOSE_PORT} -e USERPASS=${MY_PASSWORD} -v ${WORK_DIR}:/home/guest -v ${DATA_DIR}:/data $DOCKER_IMAGENAME"
echo $DOCKER_CMD
$DOCKER_CMD

printf "${SEP_STRING}URL:\t${URL}\nUsername: guest\nPassword: ${MY_PASSWORD} ${SEP_STRING}"
printf "${SEP_STRING}In RStudio:\n1. New Project: https://github.com/ibiem-2019/ibiem_2019_material.git"


printf "\n${SEP_STRING} To clean up: \n\n"
echo "rm -rf $DEMO_BASE"
echo "docker rm -f \`docker ps -aq\`"
echo "docker rmi \`docker images -aq\`"
printf "${SEP_STRING}\n"

