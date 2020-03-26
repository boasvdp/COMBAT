#!/bin/bash

set -e

if [[ -s scripts/snp-dists-alnlengths ]]
then
	echo "Script has been downloaded before and size > 0 bytes"
else
	wget -O scripts/snp-dists-alnlengths https://github.com/boasvdp/snp-dists/raw/master/snp-dists-alnlengths
	chmod u+x scripts/snp-dists-alnlengths
fi

