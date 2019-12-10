#!/bin/bash

set -e

if [[ -s maskrc-svg.py ]]
then
	echo "Script has been downloaded before and size > 0 bytes"
else
	wget https://raw.githubusercontent.com/kwongj/maskrc-svg/master/maskrc-svg.py
fi
