#ifndef MOLECULE_H
#define MOLECULE_H

#include <string>
#include <memory>

#include "helpers.h"

namespace doci {

/**
 * We will store the matrix elements of our
 * molecule in this class. It needs nice
 * orthonormal matrix elements. Currently, we supply
 * them from a HDF5 file. 
 * We use a RHF basis as starting point but you can do anything.
 */
class Molecule
{
   public:
      Molecule(std::string);

      Molecule(const Molecule &);

      Molecule(Molecule &&);

      virtual ~Molecule() = default;

      Molecule& operator=(const Molecule &);

      Molecule& operator=(Molecule &&);

      virtual double getT(int, int) const;

      virtual double getV(int, int, int, int) const;

      double get_nucl_rep() const;

      unsigned int get_n_sp() const;

      unsigned int get_n_electrons() const;

      void Print() const;

      double HF_Energy() const;

   private:
      //! store the one electron matrix elements
      std::unique_ptr<helpers::matrix> OEI;

      //! store the two electron matrix elements
      std::unique_ptr<helpers::matrix> TEI;

      //! number of electrons
      unsigned int n_electrons;

      //! nuclear repulsion (const part to add to the energy)
      double nucl_rep;

      //! the size of the single particles space (without spin)
      unsigned int n_sp;
};

}

#endif /* MOLECULE_H */

/* vim: set ts=3 sw=3 expandtab :*/
