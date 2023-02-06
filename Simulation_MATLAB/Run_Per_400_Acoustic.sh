# Author: Soroosh Sanatkhani
# Columbia University
# Created: February 2, 2023
# Last Modified: February 6, 2023
#######################################
#######################################

# Make sure to use chmod -R 777 ./ in the folder first

#!/bin/bash
cd Perpendicular_400
echo "PosN--Part 1 ..."
matlab -batch "Part1_PosN_preAcoustic" 2>&1 | tee Per_400_1_1.txt
echo "PosN--Accoustic simulation ..."
../kspaceFirstOrder-OMP -i input.h5 -s 2705 --verbose 2 -r 1% -o output.h5 -p  2>&1 | tee Per_400_1_2.txt
echo "PosN--Part 2 ..."
matlab -batch "Part2_PosN_preThermal"  2>&1 | tee Per_400_2.txt

cd ..

echo Done!