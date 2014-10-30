#ifndef DM2_H
#define DM2_H

#include <vector>
#include <memory>

#include "helpers.h"
#include "Molecule.h"
#include "Permutation.h"

/**
 * This will store an second order density matrix from a
 * DOCI wavefunction. It only stores the non-zero elements,
 * meaning: a block with dimension of the sp levels and a
 * diagonal that is four fold degenerate.
 */
class DM2
{
   public:
      DM2(unsigned int, unsigned int);

      DM2(const Molecule &);

      DM2(const DM2 &);

      DM2(DM2 &&);

      virtual ~DM2() = default;

      DM2& operator=(const DM2 &);

      DM2& operator=(DM2 &&);

      double operator()(int, int, int, int) const;

      void fill_lists(unsigned int);

      void WriteToFile(const std::string) const;

      static DM2 ReadFromFile(const std::string);

      unsigned int get_n_electrons() const;

      unsigned int get_n_sp() const;

      void Build(Permutation &, std::vector<double> &);

   private:

      //! convert single particles indices to two particles indices
      static std::unique_ptr<helpers::matrix> sp2tp;

      //! convert two particles indices to single particles indices
      static std::unique_ptr<helpers::matrix> tp2sp;

      //! the block part of the 2DM
      std::unique_ptr<helpers::matrix> block;

      //! the remaining diagonal part of the 2DM
      std::vector<double> diag;

      //! number of particles
      unsigned int n;
};

#endif /* DM2_H */

/* vim: set ts=3 sw=3 expandtab :*/
