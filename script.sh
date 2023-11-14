#!/bin/bash

cd /home/genouest/genscale/vlevallois/cqf
. /local/env/envconda.sh
conda activate ~/my_env
make test
./test 31 33
