#!/bin/bash

# This file contains a set of instructions to run the main file cuOPO3D.cu.
# Use this file to perform simulations in bulk. This means, for example, 
# systematically varying the input power, the cavity reflectivity, etc.

TIMEFORMAT='It took %0R seconds.' 
time { 

clear		# Clear screen
##  Comment next lines for safe.
rm *.dat	# This removes all .dat files in the current folder. 
rm *.txt	# This removes all .txt files in the current folder. 
rm cuOPO3D	# This removes a previuos executable file (if it exist)
# rm -r test

#####################################################################
# -----Compilation line-----
# There are different compilation modes for 3D simulations:
# 	- Define the preprocessor variable -DDIFFRACCION to include diffraction effects.
# 	- Define the preprocessor variable -DDISPERSION to include dispersion effects.
#	- Without setting any previous preprocessor variables in compilation line,
#	  code will run waveplane approximation model.


nvcc cuOPO3D.cu -DDIFFRACTION -DDISPERSION -diag-suppress 177 -gencode=arch=compute_60,code=sm_60 -gencode=arch=compute_60,code=compute_60 -O2 -lcufftw -lcufft -o cuOPO3D
# nvcc cuOPO3D.cu -DDIFFRACTION -diag-suppress 177 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_75,code=compute_75 -O2 -lcufftw -lcufft -o cuOPO3D


# There are three flags specific for CUDA compiler:
# -gencode=arch=compute_75,code=sm_75 
# -gencode=arch=compute_75,code=compute_75 : please check your GPU card architecture (Ampere, Fermi, Turing, Kepler, etc) 
#											 to set the correct number sm_XX and compute_XX.
#											 This code was tested using a Nvidia GeForce GTX 1650 card (Turing architecture). 
#
# Please visit https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#gpu-feature-list to set the proper flag according to your GPU card.
# -lcufftw -lcufft : flags needed to perform cuFFT (the CUDA Fourier transform)
#####################################################################

# The variables defined below (ARGX) will be passed as arguments to the main file 
# cuOPO3D.cu on each execution via the argv[X] instruction.

ARG1=(0.5)     	# Pump power                (ARG1)

p=$(awk 'BEGIN{for(i=0.5; i<=0.5; i=i+0.025)print i}')
w=$(awk 'BEGIN{for(i=31; i<=31; i=i+5)print i}')

# Each for-loop span over one or more values defined in the previous arguments. 
# 		for N in $x
# 		do

for WAIST in $w 
do
	FOLDERSIM="waist_${WAIST}"
	# for (( p=0; p<${#ARG1[@]}; p++ ))
	# do
	for N in $p
	do
		# N=${ARG1[$p]}
		printf "\nPower			= ${N} W\n" 

		printf "\nMaking directory...\n"
		FOLDER="sPPLT_N_${N}"
		FILE="sPPLT_N_${N}.txt"
		
		printf "Bash execution and writing output file...\n\n"
		./cuOPO3D $N $WAIST | tee -a $FILE
		
		printf "Bash finished!!\n\n" 
		mkdir $FOLDER
		mv *.dat $FOLDER"/"
		mv FILE.txt $FOLDER"/"
	done


	if [ -d "$FOLDERSIM" ]; then
		echo "Moving simulations in ${FOLDERSIM}..."
		mv sPPLT* $FOLDERSIM"/" 
		mv *.txt $FOLDERSIM"/" 
	else

		mkdir $FOLDERSIM
		echo "Creating and moving simulations in ${FOLDERSIM}..."
		mv sPPLT* $FOLDERSIM"/" 
		mv *.txt $FOLDERSIM"/" 
	fi

done

}