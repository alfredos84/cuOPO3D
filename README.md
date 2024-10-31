# cuOPO3D: GPU-based Optical Parametric Oscillator designer

is a CUDA-based toolkit for simulating optical parametric oscillators (OPOs) using the coupled-wave equations (CWEs) that well-describe the three wave mixing (TWM) processes in a second-order nonlinear media.
The code was written to be as general as possible, so that it includes both scattering and diffraction effects in a single simulation. This means the model is based on a (3+1)D physical problem (3 spatials and 1 temporal dimensions). The model includes the diffraction and the dispersion terms, and linear absorption. Intracavity element can easily be implemented if they are required. 

The simulation outcomes are the pump, signal and idler complex electric fields in the temporal, spatial and spectral domain, as required.

The user is free to use this package either in OPO or single-pass configuration (without cavity) to simulate continuous-wave, picosecond or femtosecond regimes by modifying only a few lines in the main file.

CUDA programming allows you to implement parallel computing in order to speed up calculations that typically require a considerable computational demand.

<--------------------------------------------------------------------------------------------->
## Setup and execution
For running this code is necessary to have a GPU in your computer and installed the CUDA drivers and the CUDA-TOOLKIT as well. 
To install the CUDA driver on a Linux system please visit: https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html




To run simulations using the package clone this project typing in a terminal
```
https://github.com/alfredos84/Parallel_OPO_Simulator.git
```
Once the project was cloned, the user will find a parent folder `cuOPO3D` containing two other
- `src`: contains the main file `cuOPO3D.cu`, the headers folder, and the bash file `cuOPO3D.sh` used to compile and execute the package by passing several simulations parameters.
- `src/headers`: this folder contains the headers files needed to execute the package.

### Bash file `src/cuOPO3D.sh`

The bash file is mainly used to massively perform simulations by passing the main file different parameters such as pump power, cavity detuning, etc. Before starting the user has to allow the system execute the bash file. To do that type in the terminal
```
chmod 777 cuOPO3D.sh # enable permissions for execution
```

Finally, to execute the file execute the following command line
```
./cuOPO3D.sh         # execute the files
```

In the `cuOPO.sh` file you will find the command line for the compilation:
```
nvcc cuOPO3D.cu -DDIFFRACTION -DDISPERSION -diag-suppress 177 -gencode=arch=compute_60,code=sm_60 -gencode=arch=compute_60,code=compute_60 -O2 -lcufftw -lcufft -o cuOPO3D
```
The flags `-gencode=...` are related to the GPU architecture (see below more information about this point). The flag `-O2` is the level of optimization, and can be omited or changed by the user. The flags `-lcufftw` and `-lcufft` tell the compiler to use the `CUFFT library` that performs the Fourier transform on GPU.

Finally, the execution is done using the command line in the `cuOPO3D.sh` file is
```
./cuOPO3D <ARGUMENTS_TO_PASS>
```
where `ARGUMENTS_TO_PASS` are variables externaly passed to the main file `cuOPO3D.cu`.

### Outputs

This package returns a set of `.dat` files with the signal, idler and pump electric fields, separated into real and imaginary parts. It also returns time and frequency vectors

### GPU architecture
Make sure you know your GPU architecture before compiling and running simulations. For example, pay special attention to the sm_75 flag defined in the provided `cuOPO.sh` file. That flag might not be the same for your GPU since it corresponds to a specific architecture. For instance, I tested this package using two different GPUs:
1. Nvidia Geforce MX250: architecture -> Pascal -> flag: sm_60
2. Nvidia Geforce GTX1650: architecture -> Turing -> flag: sm_75

Please check the NVIDIA documentation in https://docs.nvidia.com/cuda/


### Contact me
For any questions or queries, do not hesitate to contact the developer 
