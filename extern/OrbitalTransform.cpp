// CIFlow is a very flexible configuration interaction program
// Copyright (C) 2014 Mario Van Raemdonck <mario.vanraemdonck@UGent.be>
//
// This file is part of CIFlow.
//
// CIFlow is private software; developed by Mario Van Raemdonck
// a member of the Ghent Quantum Chemistry Group (Ghent University).
// See also : http://www.quantum.ugent.be
//
// At this moment CIFlow is not yet distributed.
// However this might change in the future in the hope that
// it will be useful to someone.
//
// For now you have to ask the main author for permission.
//
//--
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cassert>

#include "OrbitalTransform.h"
#include "Hamiltonian.h" 
#include "OptIndex.h"
#include "UnitaryMatrix.h"
#include "Lapack.h"

using std::min;
using std::max;
using CheMPS2::Hamiltonian;
using namespace simanneal;

OrbitalTransform::OrbitalTransform(const Hamiltonian& HamIn):
    index(HamIn)
{
    OptIndex index(HamIn);
    numberOfIrreps = index.getNirreps();
    int L = HamIn.getL();
    SymmInfo.setGroup(HamIn.getNGroup());
    _hamorig.reset(new CheMPS2::Hamiltonian(HamIn));

    _unitary.reset(new UnitaryMatrix(index));

    //Create the memory for the orbital transformations.
    unsigned long long maxlinsize = 0;
    for (int irrep=0; irrep< index.getNirreps(); irrep++)
    {
        unsigned int linsize_irrep = index.getNORB(irrep);
        if (linsize_irrep > maxlinsize)  
            maxlinsize  = linsize_irrep;
    }

    //Determine the blocksize for the 2-body transformation
    auto& maxBlockSize = maxlinsize;
    //Allocate 2-body rotation memory: One array is approx (maxBlockSize/273.0)^4 * 42 GiB --> [maxBlockSize=100 --> 750 MB]
    auto maxBSpower4 = maxBlockSize * maxBlockSize * maxBlockSize * maxBlockSize; //Note that 273**4 overfloats the 32 bit integer!!!
    auto sizeWorkmem1 = max( max( maxBSpower4 , 3*maxlinsize*maxlinsize ) , 1uLL * L * L ); //For (2-body tfo , updateUnitary, calcNOON)
    auto sizeWorkmem2 = max( max( maxBSpower4 , 2*maxlinsize*maxlinsize ) , L*(L + 1uLL) ); //For (2-body tfo, updateUnitary and rotate_to_active_space, rotate2DMand1DM)
    mem1.reset(new double[sizeWorkmem1]);
    mem2.reset(new double[sizeWorkmem2]);
}

void OrbitalTransform::fillHamCI(Hamiltonian& HamCI)
{
    assert(&HamCI != _hamorig.get());	
    buildOneBodyMatrixElements();	
    fillConstAndTmat(HamCI); //fill one body terms and constant part.	

    //Two-body terms --> use eightfold permutation symmetry in the irreps :-)
    for (int irrep1 = 0; irrep1<numberOfIrreps; irrep1++)
        for (int irrep2 = irrep1; irrep2<numberOfIrreps; irrep2++)
        {
            const int productSymm = SymmInfo.directProd(irrep1,irrep2);
            for (int irrep3 = irrep1; irrep3<numberOfIrreps; irrep3++)
            {
                const int irrep4 = SymmInfo.directProd(productSymm,irrep3);
                // Generated all possible combinations of allowed irreps
                if (irrep4>=irrep2)
                {
                    int linsize1 = index.getNORB(irrep1);
                    int linsize2 =index.getNORB(irrep2);
                    int linsize3 =index.getNORB(irrep3);
                    int linsize4 =index.getNORB(irrep4);

                    if ((linsize1>0) && (linsize2>0) && (linsize3>0) && (linsize4>0))
                    {
                        for (int cnt1=0; cnt1<linsize1; cnt1++)
                            for (int cnt2=0; cnt2<linsize2; cnt2++)
                                for (int cnt3=0; cnt3<linsize3; cnt3++)
                                    for (int cnt4=0; cnt4<linsize4; cnt4++)
                                        mem1[cnt1 + linsize1 * ( cnt2 + linsize2 * (cnt3 + linsize3 * cnt4) ) ]
                                            = _hamorig->getVmat(index.getNstart(irrep1) + cnt1,index.getNstart(irrep2) + cnt2, index.getNstart(irrep3) + cnt3, index.getNstart(irrep4) + cnt4 );

                        char trans = 'T';
                        char notra = 'N';
                        double alpha = 1.0;
                        double beta  = 0.0; //SET !!!

                        int rightdim = linsize2 * linsize3 * linsize4; //(ijkl) -> (ajkl)
                        double * Umx = _unitary->getBlock(irrep1);
                        dgemm_(&notra, &notra, &linsize1, &rightdim, &linsize1, &alpha, Umx, &linsize1, mem1.get(), &linsize1, &beta, mem2.get(), &linsize1);

                        int leftdim = linsize1 * linsize2 * linsize3; //(ajkl) -> (ajkd)
                        Umx = _unitary->getBlock(irrep4);
                        dgemm_(&notra, &trans, &leftdim, &linsize4, &linsize4, &alpha, mem2.get(), &leftdim, Umx, &linsize4, &beta, mem1.get(), &leftdim);

                        int jump1 = linsize1 * linsize2 * linsize3; //(ajkd) -> (ajcd)
                        int jump2 = linsize1 * linsize2 * linsize3;
                        leftdim   = linsize1 * linsize2;
                        Umx = _unitary->getBlock(irrep3);
                        for (int bla=0; bla<linsize4; bla++)
                            dgemm_(&notra, &trans, &leftdim, &linsize3, &linsize3, &alpha, mem1.get()+jump1*bla, &leftdim, Umx, &linsize3, &beta, mem2.get()+jump2*bla, &leftdim);

                        jump2    = linsize1 * linsize2;
                        jump1    = linsize1 * linsize2;
                        rightdim = linsize3 * linsize4;
                        Umx = _unitary->getBlock(irrep2);
                        for (int bla=0; bla<rightdim; bla++)
                            dgemm_(&notra, &trans, &linsize1, &linsize2, &linsize2, &alpha, mem2.get()+jump2*bla, &linsize1, Umx, &linsize2, &beta, mem1.get()+jump1*bla, &linsize1);

                        for (int cnt1=0; cnt1<linsize1; cnt1++)
                            for (int cnt2=0; cnt2<linsize2; cnt2++)
                                for (int cnt3=0; cnt3<linsize3; cnt3++)
                                    for (int cnt4=0; cnt4<linsize4; cnt4++)
                                        HamCI.setVmat(index.getNstart(irrep1) + cnt1,index.getNstart(irrep2) + cnt2, index.getNstart(irrep3) + cnt3,index.getNstart(irrep4) + cnt4, mem1[cnt1 + linsize1 * ( cnt2 + linsize2 * (cnt3 + linsize3 * cnt4) ) ] );

                    } //end if the problem has orbitals from all 4 selected irreps
                } // end if irrep 4 >= irrep2
            }// end run irrep3
        } // end run irrep2
}

/**
 * This method fills the OneBodyMatrixElements array with
 * the 1-body matrix elements so we can rotate them
 */
void OrbitalTransform::buildOneBodyMatrixElements()
{
    if(!QmatrixWork.size() || !OneBodyMatrixElements.size())
    {
        QmatrixWork.resize(numberOfIrreps);
        OneBodyMatrixElements.resize(numberOfIrreps);

        for (int irrep=0; irrep<numberOfIrreps; irrep++)
        {
            const int size = index.getNORB(irrep) * index.getNORB(irrep);
            QmatrixWork[irrep].reset(new double[size]);
            OneBodyMatrixElements[irrep].reset(new double[size]);
        }	
    }

    //#pragma omp parallel for schedule(dynamic)
    for (int irrep=0; irrep<numberOfIrreps; irrep++)
    {
        const int linsize = index.getNORB(irrep);
        for (int row=0; row<linsize; row++)
        {
            const int HamIndex1 = index.getNstart(irrep) + row;
            //set the diagonal elements.
            OneBodyMatrixElements[irrep][row*(1+linsize)] = _hamorig->getTmat(HamIndex1,HamIndex1);
            for (int col=row+1; col<linsize; col++)
            {
                const int HamIndex2 = index.getNstart(irrep) + col;
                // remember blas and lapack are in fortran and there it is convenient to run first over columnsand then over rows.
                OneBodyMatrixElements[irrep][row + linsize * col] = _hamorig->getTmat(HamIndex1,HamIndex2); 
                OneBodyMatrixElements[irrep][col + linsize * row] = OneBodyMatrixElements[irrep][row + linsize * col];
            }
        }
    }
    rotate_old_to_new(OneBodyMatrixElements.data());
}

/**
 * Calculates the rotation of the 1-body matrix
 * T' = Q * T * Q^T
 * It uses the OneBodyMatrixElements array as start and endpoint
 * point and uses the Qmatrixwork as temp storage
 */
void OrbitalTransform::rotate_old_to_new(std::unique_ptr<double []> * matrix)
{
    //#pragma omp parallel for schedule(dynamic)
    for (int irrep=0; irrep<numberOfIrreps; irrep++)
    {
        int linsize  = get_norb(irrep);
        if(linsize > 1)
        {
            double * Umx = _unitary->getBlock(irrep);
            double alpha = 1.0;
            double beta  = 0.0;
            char trans   = 'T';
            char notrans = 'N';
            // Qmatrixwork = Umx * OneBodyMatrixElements
            dgemm_(&notrans,&notrans,&linsize,&linsize,&linsize,&alpha,Umx,&linsize,matrix[irrep].get(),&linsize,&beta,QmatrixWork[irrep].get(),&linsize);
            // OneBodyMatrixElements = Qmatrixwork * Umx^T
            dgemm_(&notrans,&trans,  &linsize,&linsize,&linsize,&alpha,QmatrixWork[irrep].get(),&linsize,Umx,&linsize,&beta,matrix[irrep].get(),&linsize);
        }
    }
}

double OrbitalTransform::TmatRotated(const int index1, const int index2) const
{
    if(!OneBodyMatrixElements.size())
    {
        std::cerr << "First build the OneBodyMatrixElements!" << std::endl;
        exit(EXIT_FAILURE);
    }

    const int irrep1 = _hamorig->getOrbitalIrrep(index1);
    const int irrep2 = _hamorig->getOrbitalIrrep(index2);

    if (irrep1 != irrep2)
        //From now on: both irreps are the same.
        return 0.0;

    int shift = index.getNstart(irrep1);

    return OneBodyMatrixElements[irrep1][index1-shift + index.getNORB(irrep1) * (index2-shift)];
}

void OrbitalTransform::fillConstAndTmat(Hamiltonian& Ham) const
{
    //Constant part of the energy
    double value = _hamorig->getEconst();
    Ham.setEconst(value);

    //One-body terms: diagonal in the irreps
    for (int irrep=0; irrep<numberOfIrreps; irrep++)
    {
        const int linsize = index.getNORB(irrep);
        for (int cnt1=0; cnt1<linsize; cnt1++)
            for (int cnt2=cnt1; cnt2<linsize; cnt2++)
            {
                int shift = index.getNstart(irrep);
                Ham.setTmat( cnt1 + shift , cnt2 + shift , TmatRotated(cnt1 + shift , cnt2 + shift) );
            }
    }
}

double OrbitalTransform::get_difference_orig(Hamiltonian& hamin) const
{
    assert(0);
    return 0;
//    return _hamorig->get_difference(&hamin);
}

void OrbitalTransform::CheckDeviationFromUnitary() const
{
    _unitary->CheckDeviationFromUnitary(mem2.get());
}

void OrbitalTransform::set_unitary(UnitaryMatrix& unit)
{
    _unitary.reset(new UnitaryMatrix(unit));
}

double OrbitalTransform::get_norb(int irrep) const
{
    return index.getNORB(irrep);
}

void OrbitalTransform::update_unitary(double * change)
{
    _unitary->updateUnitary(mem1.get(),mem2.get(),change,1);//First 1, sets multiply.
}

/* vim: set ts=4 sw=4 expandtab :*/
