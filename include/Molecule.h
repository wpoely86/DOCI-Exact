#ifndef MOLECULE_H
#define MOLECULE_H

#include <string>
#include <memory>

#include "helpers.h"

namespace doci {

/**
 * This class is the inferface for all
 * integrals classes. It's pure virtual.
 */
class Molecule
{
   public:

      virtual ~Molecule() = default;

      virtual Molecule* clone() const = 0;

      virtual double getT(int, int) const = 0;

      virtual double getV(int, int, int, int) const = 0;

      virtual double get_nucl_rep() const;

      virtual unsigned int get_n_sp() const;

      virtual unsigned int get_n_electrons() const;

   protected:

      //! number of electrons
      unsigned int n_electrons;

      //! nuclear repulsion (const part to add to the energy)
      double nucl_rep;

      //! the size of the single particles space (without spin)
      unsigned int n_sp;
};


/**
 * We will store the matrix elements generated
 * by our PSI4 plugin in this this class. We don't use any
 * symmetry (aka, C1 symmetry). It needs nice orthonormal
 * matrix elements. Currently, we supply them from a HDF5 file. 
 * We use a RHF basis as starting point but you can do anything.
 */
class PSI_C1_Molecule: public Molecule
{
   public:
      PSI_C1_Molecule(std::string);

      PSI_C1_Molecule(const PSI_C1_Molecule &);

      PSI_C1_Molecule(PSI_C1_Molecule &&);

      PSI_C1_Molecule& operator=(const PSI_C1_Molecule &);

      PSI_C1_Molecule& operator=(PSI_C1_Molecule &&);

      PSI_C1_Molecule* clone() const;

      double getT(int, int) const;

      double getV(int, int, int, int) const;

      void Print() const;

      double HF_Energy() const;

   private:
      //! store the one electron matrix elements
      std::unique_ptr<helpers::matrix> OEI;

      //! store the two electron matrix elements
      std::unique_ptr<helpers::matrix> TEI;
};

}

#endif /* MOLECULE_H */

/* vim: set ts=3 sw=3 expandtab :*/
