/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013, 2014 Sebastian Wouters

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along
   with this program; if not, write to the Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef OPTIONS_CHEMPS2_H
#define OPTIONS_CHEMPS2_H

#include <stdlib.h>
#include <string>

using std::string;

namespace CheMPS2{

   const int    DMRGSCF_maxIterations         = 100;
   const double DMRGSCF_gradientNormThreshold = 1e-6;
   const bool   DMRGSCF_storeUnitary          = true;
   const string DMRGSCF_unitaryStorageName    = "CheMPS2_CASSCF.h5";
   const int    DMRGSCF_maxlinsizeCutoff      = 100;
   const bool   DMRGSCF_debugPrint            = false;
   const bool   DMRGSCF_stateAveraged         = true;
   
   const int    DMRGSCF_whichActiveSpace      = 0;
   const bool   DMRGSCF_dumpCorrelations      = false;
   const bool   DMRGSCF_startLocRandom        = false;
   
   const bool   DMRGSCF_doDIIS                = false;
   const double DMRGSCF_DIISgradientBranch    = 1e-2;
   const int    DMRGSCF_numDIISvecs           = 7;
   const bool   DMRGSCF_storeDIIS             = true;
   const string DMRGSCF_DIISstorageName       = "CheMPS2_DIIS.h5";

   const string TMPpath                       = "/tmp";
   const bool   DMRG_printDiscardedWeight     = false;
   const bool   DMRG_storeRenormOptrOnDisk    = true;
   const bool   DMRG_storeMpsOnDisk           = false;
   const string DMRG_MPS_storage_prefix       = "CheMPS2_MPS";
   const string DMRG_OPERATOR_storage_prefix  = "CheMPS2_Operators_";
   
   const bool   HAMILTONIAN_debugPrint        = false;
   const string HAMILTONIAN_TmatStorageName   = "CheMPS2_Ham_Tmat.h5";
   const string HAMILTONIAN_VmatStorageName   = "CheMPS2_Ham_Vmat.h5";
   const string HAMILTONIAN_ParentStorageName = "CheMPS2_Ham_parent.h5";
   
   const string TWODM_2DM_A_storagename       = "CheMPS2_2DM-A.h5";
   const string TWODM_2DM_B_storagename       = "CheMPS2_2DM-B.h5";
   
   const bool   HEFF_debugPrint               = true;
   const int    HEFF_DAVIDSON_NUM_VEC         = 32;
   const int    HEFF_DAVIDSON_NUM_VEC_KEEP    = 3;
   const double HEFF_DAVIDSON_PRECOND_CUTOFF  = 1e-12;
   const double HEFF_DAVIDSON_RTOL_BASE       = 1e-10;
   
   const bool   SYBK_debugPrint               = false;
   const int    SYBK_dimensionCutoff          = 262144;
   
   const double TENSORT_orthoComparison       = 1e-13;
   
   const bool   CORRELATIONS_debugPrint       = false;
   const double CORRELATIONS_discardEig       = 1e-100;
   
   const double EDMISTONRUED_gradThreshold    = 1e-8;
   const int    EDMISTONRUED_maxIter          = 1000;
   const int    EDMISTONRUED_maxIterBackTfo   = 15;

   const bool   Orbopt_debugPrint             = false;
}

#endif

