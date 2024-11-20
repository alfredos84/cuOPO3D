#include "DtypesConsts.cuh"
#include "Common.cuh"
#include "Files.cuh"
#include "Operators.cuh"
#include "Crystals.cuh"

#ifdef DIFFRACTION
#include "EFields3D.cuh"
#include "Cavity3D.cuh"
#include "Solver3D.cuh"
#else
#include "EFields1D.cuh"
#include "Cavity1D.cuh"
#include "Solver1D.cuh"
#endif