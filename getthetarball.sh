#!/bin/bash

scp -i ~/.ssh/id_rsa neilmacl@vortex.ccr.buffalo.edu:data.tar ~/Documents/early-warning-spatial

cd ./data/
mv *.rds ./old/
cd ../

mv ./data.tar ./data/
cd ./data/
tar -xvf data.tar
ls -l | head -n 10
