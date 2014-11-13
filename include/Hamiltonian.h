#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <vector>
#include <memory>

#include "Permutation.h"
#include "Molecule.h"
#include "SparseMatrix_CRS.h"

namespace doci {

/**
 * DOCIHamiltonian will store the actual hamiltonian in a sparse format
 * It needs a Permutation object for the basisset and a Molecule object
 * for the matrix elements
 */
class DOCIHamiltonian
{
   public:
      DOCIHamiltonian(const Permutation &,const Molecule &);

      DOCIHamiltonian(const Molecule &);

      DOCIHamiltonian(const DOCIHamiltonian &);

      DOCIHamiltonian(DOCIHamiltonian &&) = default;

      virtual ~DOCIHamiltonian() = default;

      DOCIHamiltonian& operator=(const DOCIHamiltonian &);

      DOCIHamiltonian& operator=(DOCIHamiltonian &&) = default;

      Molecule const & getMolecule() const;

      Permutation const & getPermutation() const;

      unsigned int getdim() const;

      void Build();

      std::pair< std::vector<double>,helpers::matrix > DiagonalizeFull() const;

      std::pair< double,std::vector<double> > Diagonalize() const;

      static unsigned int CountBits(mybitset);

      static int CalcSign(unsigned int i,unsigned int j, const mybitset a);

   private:

      void Build_iter(Permutation &, helpers::SparseMatrix_CRS &,unsigned int, unsigned int, Molecule &);

      std::unique_ptr<Permutation> permutations;

      std::unique_ptr<Molecule> molecule;

      std::unique_ptr<helpers::SparseMatrix_CRS> mat;
};

}

#endif /* HAMILTONIAN_H */

/* vim: set ts=3 sw=3 expandtab :*/
