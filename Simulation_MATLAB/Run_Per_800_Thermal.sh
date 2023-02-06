# Author: Soroosh Sanatkhani
# Columbia University
# Created: February 2, 2023
# Last Modified: February 6, 2023
#######################################
#######################################

# Make sure to use chmod -R 777 ./ in the folder first

#!/bin/bash
cd Perpendicular_800
echo "PosN--Part 3 - Thermal simulation ..."
matlab -batch "Part3_PosN_Thermal" 2>&1 | tee Per_800_3.txt

cd ..

echo Done!