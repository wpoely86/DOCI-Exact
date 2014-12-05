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
#ifndef _UnitaryMatrix_H
#define _UnitaryMatrix_H

#include <memory>
#include <vector>
#include "OptIndex.h" //For optindex (should put this in a separate file)

namespace simanneal {

/** unitary class.
    The UnitaryMatrix class is a storage and manipulation class for the unitary matrix. This matrix is blockdiagonal in the irreducible representations, and is formed by stepwise multiplying in new unitary rotations due to the Newton-Raphson algorithm.
*/

class UnitaryMatrix
{
    public:
      
        //! Constructor
        /** \param _hamindexIn */
        UnitaryMatrix(OptIndex& index);

     	//!Copy constructor
        UnitaryMatrix(const UnitaryMatrix & unit);

        UnitaryMatrix(UnitaryMatrix && unit);

        //! Destructor
        virtual ~UnitaryMatrix() = default;

        UnitaryMatrix& operator=(const UnitaryMatrix &unit);

        UnitaryMatrix& operator=(UnitaryMatrix &&unit);
        
        //! Get the number of variables in the x-parametrization of the unitary update
        /** \return The number of unique variables in the x-matrix */
        unsigned int getNumVariablesX() const;
        
        //! Get the first Hamiltonian index corresponding to linearindex
        /** \param linearindex The linear index of the x-parametrization
            \return The first Hamiltonian index corresponding to linearindex */
        int getFirstIndex(const int linearindex) const;
        
        //! Get the second Hamiltonian index corresponding to linearindex
        /** \param linearindex The linear index of the x-parametrization
            \return The second Hamiltonian index corresponding to linearindex */
        int getSecondIndex(const int linearindex) const;
        
        //! Get the unitary rotation for block irrep
        /** \param irrep The irreducible representation
            \return Pointer to the desired unitary block */
        double * getBlock(const int irrep) const;
        
        //! Copy the x solution back (NR, augmented NR, ...)
        /** \param vector The x-solution */
        void copyXsolutionBack(double * vector);
        
        //! Update the unitary transformation based on the new vector and the previous unitary
        /** \param workmem1 Work memory
            \param workmem2 Work memory */
        void updateUnitary(double * workmem1, double * workmem2, double * vector , const bool multiply);
        
        //! Rotate the unitary matrix to the NO eigenbasis
        /** \param eigenvecs The NO eigenbasis
            \param work Work memory */
        void rotate_active_space_vectors(double * eigenvecs, double * work);
        
        //! Calculate the two-norm of U^T*U - I
        /** \param work Work memory */
        void CheckDeviationFromUnitary(double * work) const;
        
        //! Save the unitary to disk
        void saveU(std::string savename) const;
        
        //! Load the unitary from disk
        void loadU(std::string loadname);
        
        //! Delete the stored unitary (on disk)
        void deleteStoredUnitary(std::string name) const;

        //! Implements a simple jacobi rotation of the i, and j th orbital.
        void jacobi_rotation(int irrep, int i, int j, double angle);

        //! Resets the unitary matrix
        void reset_unitary();

        void print_unitary() const;

        void build_skew_symm_x(const int irrep, double * xblock , const double * Xelem) const;

        void sendreceive(int);

        int get_Nirrep() const;

    private:
      
        //Externally created and destroyed index handler
        std::unique_ptr<OptIndex> _index;
      
        //Number of variables in the x-matrix
        unsigned int x_linearlength;
        
        //The unitary matrix (e^x * previous unitary): unitary[irrep][row + size_irrep * col]
        std::vector< std::unique_ptr<double []> > unitary;
        
        // Find the linear index corresponding to p and q
        /** \param p_index The first Hamiltonian index
            \param q_index The second Hamiltonian index
            \return The linear index corresponding to (p,q). If no index is found -1 is returned. */
         
   };

}
#endif

/* vim: set ts=4 sw=4 expandtab :*/
