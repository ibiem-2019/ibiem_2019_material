
SEP_STRING="\n--------------------------------------------------\n"

REMOTE_OR_LOCAL=${1:-"LOCAL"}
MY_PASSWORD=${2:-"BadPa55word"}


# echo ${1}
# echo "---${REMOTE_OR_LOCAL}--"
# echo "---${1}--"


HOST_PORT=8787

if [ $REMOTE_OR_LOCAL == "REMOTE" ]; then
    EXPOSE_PORT="${HOST_PORT}:8787"
    URL="http://`hostname -A | cut -f1 -d' '`:${HOST_PORT}/"
else
    EXPOSE_PORT="127.0.0.1\:${HOST_PORT}\:8787"
    URL="http://localhost:${HOST_PORT}/"
fi




echo "expose argument: ${EXPOSE_PORT}"

printf "${SEP_STRING}URL:\t${URL}\nUsername: guest\nPassword: ${MY_PASSWORD} ${SEP_STRING}"
