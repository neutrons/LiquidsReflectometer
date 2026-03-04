#!/usr/bin/bash

# import library to do the real work
# shellcheck disable=SC1091
#. "$(dirname "$(realpath "$0")")/nsd-app-wrap.sh" \
#  || . /bin/nsd-app-wrap.sh
. /bin/nsd-app-wrap.sh


# put together arguments - env, application, argv
args=("lr_reduction_dev" "python /SNS/REF_L/shared/launcher/launcher.py" "$@")
# launch the tool
# BUG: needs to be in the /SNS/REF_L/shared/launcher/ when starting
( cd /SNS/REF_L/shared/launcher && pixi_launch "${args[@]}" )
