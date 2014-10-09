#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "Permutation.h"
#include "Molecule.h"
#include "SparseMatrix_CCS.h"

/**
 * Hamiltonian will store the actual hamiltonian in a sparse format
 * It needs a Permutation object for the basisset and a Molecule object
 * for the matrix elements
 */
class Hamiltonian
{
   public:
      Hamiltonian(Permutation &, Molecule &);

      Hamiltonian(Molecule &);


      Molecule const & getMolecule() const;

      Permutation const & getPermutation() const;

      void Build();

   private:

      Permutation perm;

      Molecule molecule;

};

#endif /* HAMILTONIAN_H */

/* vim: set ts=3 sw=3 expandtab :*/
