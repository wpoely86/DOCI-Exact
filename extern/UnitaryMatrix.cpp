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
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cstring>
//#include <mpi.h>

#include "MyHDF5.h"
#include "Lapack.h"
#include "UnitaryMatrix.h"
#include "Options.h"

using std::string;
using std::ifstream;
using std::cerr;
using std::cout;
using std::endl;
using CheMPS2::Orbopt_debugPrint;
using namespace simanneal;

UnitaryMatrix::UnitaryMatrix(OptIndex& index)
{
    _index.reset(new OptIndex(index));

    //Allocate the unitary
    unitary.resize(_index->getNirreps());
    x_linearlength = 0;
    for (int irrep=0; irrep<_index->getNirreps(); irrep++)
    {
        const int linsize = _index->getNORB(irrep);
        const int size = linsize * linsize;
        x_linearlength += size;
        unitary[irrep].reset(new double[size]);
        memset(unitary[irrep].get(), 0, sizeof(double)*size);
        for (int cnt=0; cnt<linsize; cnt++)
            unitary[irrep][cnt*(1+linsize)] = 1.0;
    }
}

UnitaryMatrix::UnitaryMatrix(const UnitaryMatrix & unit)
{
    _index.reset(new OptIndex(*unit._index));
    //Copy the unitary
    unitary.resize(_index->getNirreps());
    x_linearlength = unit.x_linearlength;

    for (int irrep=0; irrep<_index->getNirreps(); irrep++)
    {
        const int linsize = _index->getNORB(irrep);
        const int size = linsize * linsize;
        unitary[irrep].reset(new double[size]);
        for (int index=0; index<size; index++)
            unitary[irrep][index] = unit.getBlock(irrep)[index];
    }
}

UnitaryMatrix::UnitaryMatrix(UnitaryMatrix && unit)
{
    _index = std::move(unit._index);
    unitary = std::move(unit.unitary);
}

UnitaryMatrix& UnitaryMatrix::operator=(const UnitaryMatrix &unit)
{
    _index.reset(new OptIndex(*unit._index));
    //Copy the unitary
    unitary.resize(_index->getNirreps());
    x_linearlength = unit.x_linearlength;

    for (int irrep=0; irrep<_index->getNirreps(); irrep++)
    {
        const int linsize = _index->getNORB(irrep);
        const int size = linsize * linsize;
        unitary[irrep].reset(new double[size]);
        for (int index=0; index<size; index++)
            unitary[irrep][index] = unit.getBlock(irrep)[index];
    }

    return *this;
}

UnitaryMatrix& UnitaryMatrix::operator=(UnitaryMatrix &&unit)
{
    _index = std::move(unit._index);
    unitary = std::move(unit.unitary);

    return *this;
}

void UnitaryMatrix::jacobi_rotation(int irrep , int i , int j , double angle)
{
    //Simple implementation of jacobi rotation. 
    i -= _index->getNstart(irrep);
    j -= _index->getNstart(irrep);
    int linsize = _index->getNORB(irrep);
    double  work[2*linsize];

    //cout << i << j<< linsize << cos(angle) << sin(angle);
    for(int l = 0 ; l < linsize ; l++)
    {
        work[l] = unitary[irrep][l+i*linsize]*cos(angle) + unitary[irrep][l+j*linsize] * sin(angle);
        work[l+linsize] = unitary[irrep][l+j*linsize]*cos(angle) + unitary[irrep][l+i*linsize] * sin(angle)*-1.;
    }

    for(int l = 0 ; l < _index->getNORB(irrep) ; l++)
    {
        unitary[irrep][l+i*linsize] = work[l];
        unitary[irrep][l+j*linsize] = work[l+linsize];
    }
}

void UnitaryMatrix::reset_unitary()
{
   for (int irrep=0; irrep<_index->getNirreps(); irrep++)
   {
      const int linsize = _index->getNORB(irrep);
      const int size = linsize * linsize;
      memset(unitary[irrep].get(), 0, sizeof(double)*size);
      for (int cnt=0; cnt<linsize; cnt++)
          unitary[irrep][cnt*(1+linsize)] = 1.0;
   }
}

unsigned int UnitaryMatrix::getNumVariablesX() const{ return x_linearlength; }

double * UnitaryMatrix::getBlock(const int irrep) const { return unitary[irrep].get(); }

void UnitaryMatrix::updateUnitary(double * temp1, double * temp2, double * vector, const bool multiply)
{
    //Per irrep
    for (int irrep=0; irrep<_index->getNirreps(); irrep++)
    {
        int linsize = _index->getNORB(irrep);
        int size = linsize * linsize;

        //linsize is op z'n minst 2 dus temp1, temp1+size, temp1+2*size,temp1+3*size zijn zeker ok
        if (linsize>1)
        {         
            //Construct the anti-symmetric x-matrix
            double * xblock    = temp1;
            double * Bmat      = temp1 +   size;   //linsize*linsize
            double * work1     = temp1 + 2*size;   //linsize*linsize
            double * work2     = temp1 + 3*size;   //linsize*linsize

            double * workLARGE = temp2;  //4*size
            int     lworkLARGE = 4*size; //4*size = 4*linsize*linsize > 3*linsize-1

            //Construct the antisymmetric x-matrix
            build_skew_symm_x(irrep, xblock , vector);

            //Bmat <= xblock * xblock
            char notr = 'N';
            double alpha = 1.0;
            double beta = 0.0; //SET !!!
            dgemm_(&notr,&notr,&linsize,&linsize,&linsize,&alpha,xblock,&linsize,xblock,&linsize,&beta,Bmat,&linsize); 

            //Calculate its eigenvalues and eigenvectors
            //Bmat * work1 * Bmat^T <= xblock * xblock
            char uplo = 'U';
            char jobz = 'V';
            int info;
            dsyev_(&jobz, &uplo, &linsize, Bmat, &linsize, work1, workLARGE, &lworkLARGE, &info); 

            if (info != 0 )
            {
                cerr << "A problem occured in the diagonalisation of the anti-symmetric matrix squared to generated the unitary -> exp(X).";
                exit(EXIT_FAILURE);	 
            }

            //work2 <= Bmat^T * xblock * Bmat
            dgemm_(&notr,&notr,&linsize,&linsize,&linsize,&alpha,xblock,&linsize,Bmat,&linsize,&beta,work1,&linsize);
            char trans = 'T';
            dgemm_(&trans,&notr,&linsize,&linsize,&linsize,&alpha,Bmat,&linsize,work1,&linsize,&beta,work2,&linsize); //work2 = Bmat^T * xblock * Bmat

            if (Orbopt_debugPrint)
            {
                cout << "   UnitaryMatrix::updateUnitary : Lambdas of irrep block " << irrep << " : " << endl;
                for (int cnt=0; cnt<linsize/2; cnt++)
                {
                    cout << "      block = [ " << work2[2*cnt   + linsize*2*cnt] << " , " << work2[2*cnt   + linsize*(2*cnt+1)] << " ] " << endl;
                    cout << "              [ " << work2[2*cnt+1 + linsize*2*cnt] << " , " << work2[2*cnt+1 + linsize*(2*cnt+1)] << " ] " << endl;
                }
            }

            //work1 <= values of the antisymmetric 2x2 blocks
            for (int cnt=0; cnt<linsize/2; cnt++)
                work1[cnt] = 0.5*( work2[2*cnt + linsize*(2*cnt+1)] - work2[2*cnt+1 + linsize*(2*cnt)] );

            if (Orbopt_debugPrint)
            {
                for (int cnt=0; cnt<linsize/2; cnt++)
                {
                    work2[2*cnt + linsize*(2*cnt+1)] -= work1[cnt];
                    work2[2*cnt+1 + linsize*(2*cnt)] += work1[cnt];
                }

                int inc = 1;
                double RMSdeviation = ddot_(&size, work2, &inc, work2, &inc);
                RMSdeviation = sqrt(RMSdeviation);

                cout << "   UnitaryMatrix::updateUnitary : RMSdeviation of irrep block " << irrep << " (should be 0.0) = " << RMSdeviation << endl;
            }

            //Calculate exp(x)
            //work2 <= exp(Bmat^T * xblock * Bmat)
            memset(work2, 0, sizeof(double)*size);
            for (int cnt=0; cnt<linsize/2; cnt++)
            {
                double cosine = cos(work1[cnt]);
                double sine = sin(work1[cnt]);
                work2[2*cnt   + linsize*(2*cnt  )] = cosine;
                work2[2*cnt+1 + linsize*(2*cnt+1)] = cosine;
                work2[2*cnt   + linsize*(2*cnt+1)] = sine;
                work2[2*cnt+1 + linsize*(2*cnt  )] = - sine;
            }

            //Handles the case when there is an odd number of orbitals.
            for (int cnt=2*(linsize/2); cnt<linsize; cnt++)
                work2[cnt*(linsize + 1)] = 1.0;

            //work2 <= Bmat * exp(Bmat^T * xblock * Bmat) * Bmat^T = exp(xblock)
            dgemm_(&notr,&notr,&linsize,&linsize,&linsize,&alpha,Bmat,&linsize,work2,&linsize,&beta,work1,&linsize);
            dgemm_(&notr,&trans,&linsize,&linsize,&linsize,&alpha,work1,&linsize,Bmat,&linsize,&beta,work2,&linsize); //work2 = exp(xblock)

            //U <-- exp(x) * U
            int inc = 1;
            if (multiply)
            {
                //U <-- exp(x) * U
                dgemm_(&notr,&notr,&linsize,&linsize,&linsize,&alpha,work2,&linsize,unitary[irrep].get(),&linsize,&beta,work1,&linsize);
                dcopy_(&size, work1, &inc, unitary[irrep].get(), &inc);
            }
            else  //U <-- exp(x)
                dcopy_(&size, work2, &inc, unitary[irrep].get(), &inc);
        }
    }

    if (Orbopt_debugPrint){ CheckDeviationFromUnitary(temp2); }
}

void UnitaryMatrix::build_skew_symm_x(const int irrep, double * result , const double * Xelem) const
{
    //should have size: linsize*linsize.	
    //Remark we build the antisymmetrix X matrix within an irrep.
    const int linsize = _index->getNORB(irrep);

    int jump = 0;
    for (int cnt=0; cnt<irrep; cnt++)
    {
        int linsizeCNT = _index->getNORB(cnt);
        jump += linsizeCNT * (linsizeCNT-1) / 2;
    }

    for (int row=0; row<linsize; row++)
    {
        result[ row + linsize * row ] = 0;

        for (int col=row+1; col<linsize; col++)
        {
            result[ row + linsize * col ] =   Xelem[ jump + row + col*(col-1)/2 ];
            result[ col + linsize * row ] = - Xelem[ jump + row + col*(col-1)/2 ];
        }
    }
}

void UnitaryMatrix::rotate_active_space_vectors(double * eigenvecs, double * work)
{
   int passed = 0;
   int nOrbDMRG = _index->getL();
   for (int irrep=0; irrep<_index->getNirreps(); irrep++){

      const int NDMRG = _index->getNORB(irrep);
      if (NDMRG > 1){

         int rotationlinsize = NDMRG;
         int blocklinsize = _index->getNORB(irrep);

         double * temp1 = work;
         double * temp2 = work + rotationlinsize*blocklinsize;
         double * BlockEigen = eigenvecs + passed * (nOrbDMRG + 1); //after rotating the irrep before point to next eigenvector irrep block
      
         for (int row = 0; row<rotationlinsize; row++){
            for (int col = 0; col<blocklinsize; col++){
               temp1[row + rotationlinsize*col] = unitary[irrep][ row + blocklinsize * col ];
            }
         }

         char tran = 'T';
         char notr = 'N';
         double alpha = 1.0;
         double beta = 0.0;
         dgemm_(&tran,&notr,&rotationlinsize,&blocklinsize,&rotationlinsize,&alpha,BlockEigen,&nOrbDMRG,temp1,&rotationlinsize,&beta,temp2,&rotationlinsize);

         for (int row = 0; row<rotationlinsize; row++){
            for (int col = 0; col<blocklinsize; col++){
               unitary[irrep][  row + blocklinsize * col ] = temp2[row + rotationlinsize*col];
            }
         }

      }

      passed += NDMRG;

   }
   
   if (Orbopt_debugPrint){ CheckDeviationFromUnitary(work); }
}

void UnitaryMatrix::CheckDeviationFromUnitary(double * work) const
{
    char tran = 'T';
    char notr = 'N';
    double alpha = 1.0;
    double beta = 0.0;

    for (int irrep=0; irrep<_index->getNirreps(); irrep++)
    {
        int linsize = _index->getNORB(irrep);

        // linsize > 0 (if the system has orbitals within the irrep)
        if(linsize > 0)
        {
            dgemm_(&tran,&notr,&linsize,&linsize,&linsize,&alpha,unitary[irrep].get(),&linsize,unitary[irrep].get(),&linsize,&beta,work,&linsize);

            int inc = 1;
            int n = linsize*linsize;
            double norm = ddot_(&n, work, &inc, work, &inc);

            // we don't have to calculate the sqrt (expensive!) to see if we're still OK.
            norm = sqrt(norm);

            cout << "Two-norm of unitary[" << irrep << "]^(dagger) * unitary[" << irrep << "] - I = " << norm << endl;
            if((norm-1) > 1e-13)
            {
                cerr << "WARNING: we reseted the unitary because we lost unitarity." << endl;
                exit(EXIT_FAILURE);
            }
        }
    }
}

void UnitaryMatrix::saveU(string savename) const
{
    hid_t file_id = H5Fcreate(savename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t group_id = H5Gcreate(file_id, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    for (int irrep=0; irrep<_index->getNirreps(); irrep++)
    {
        int norb = _index->getNORB(irrep); 

        if(norb > 0)
        {
            std::stringstream irrepname;
            irrepname << "irrep_" << irrep;

            hsize_t dimarray      = norb * norb;
            hid_t dataspace_id    = H5Screate_simple(1, &dimarray, NULL);
            hid_t dataset_id      = H5Dcreate(group_id, irrepname.str().c_str(), H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, unitary[irrep].get());

            H5Dclose(dataset_id);
            H5Sclose(dataspace_id);
        }
    }

    H5Gclose(group_id);
    H5Fclose(file_id);
}

void UnitaryMatrix::loadU(string unitname)
{
   hid_t file_id = H5Fopen(unitname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   hid_t group_id = H5Gopen(file_id, "/Data", H5P_DEFAULT);
       
   for (int irrep=0; irrep<_index->getNirreps(); irrep++)
   {
      std::stringstream irrepname;
      irrepname << "irrep_" << irrep;

      hid_t dataset_id = H5Dopen(group_id, irrepname.str().c_str(), H5P_DEFAULT);
      H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, unitary[irrep].get());
         
      H5Dclose(dataset_id);
   }

   H5Gclose(group_id);
   H5Fclose(file_id);
}

void UnitaryMatrix::deleteStoredUnitary(std::string name) const
{
   std::stringstream temp;
   temp << "rm " << name;
   int info = system(temp.str().c_str());
   cout << "Info on CASSCF::Unitary rm call to system: " << info << endl;
}

void UnitaryMatrix::print_unitary() const
{
    ///cout << std::setprecision(5);
    for( int irrep = 0 ; irrep < _index->getNirreps() ; irrep ++)
    {
        cout << "#irrep : " << irrep << endl;
        int linsize = _index->getNORB(irrep);
        for( int j = 0 ; j < linsize ; j ++)
        {
            for( int l = 0 ; l < linsize ; l ++)
                cout << unitary[irrep][j + linsize* l] << " ";

            cout << endl;
        }
        cout << endl;
    }
}

/**
 * Spread the UnitaryMatrix from rank orig to all
 * other UnitaryMatrices in the world
 */
void UnitaryMatrix::sendreceive(int orig)
{
    for (int irrep=0; irrep<_index->getNirreps(); irrep++)
    {
        const int linsize = _index->getNORB(irrep);
        const int size = linsize * linsize;
//        MPI_Bcast(unitary[irrep].get(), size, MPI_DOUBLE, orig, MPI_COMM_WORLD);
    }
}


int UnitaryMatrix::get_Nirrep() const
{
    return _index->getNirreps();
}

/* vim: set ts=4 sw=4 expandtab :*/
