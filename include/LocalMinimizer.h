#ifndef LOCALMINIMIZER_H
#define LOCALMINIMIZER_H

#include <memory>
#include <vector>
#include <tuple>
#include <random>

#include "DOCIHamtilonian.h"
#include "SymMolecule.h"
#include "DM2.h"

#include "UnitaryMatrix.h"
#include "OrbitalTransform.h"


namespace doci {

class LocalMinimizer
{
   public:
      LocalMinimizer(const doci::Sym_Molecule &);

      LocalMinimizer(doci::Sym_Molecule &&);

      virtual ~LocalMinimizer();

      void Minimize(bool dist_choice=false);

      double get_energy() const;

      double calc_new_energy();

      double calc_new_energy(const doci::Sym_Molecule &);

      void calc_energy();

      simanneal::UnitaryMatrix& get_Optimal_Unitary();

      doci::Sym_Molecule& getHam() const;

      simanneal::OrbitalTransform& getOrbitalTf() const;

      std::vector<std::tuple<int,int,double,double>> scan_orbitals();

      double get_conv_crit() const;

      void set_conv_crit(double);

      void set_conv_steps(int);

      int choose_orbitalpair(std::vector<std::tuple<int,int,double,double>> &);

      const doci::DM2& get_DM2() const;

   private:

      //! criteria for convergence of the minimizer
      double conv_crit;

      double energy;

      //! number of steps in convergence area
      int conv_steps;

      std::unique_ptr<doci::DOCIHamiltonian> method;

      std::unique_ptr<doci::DM2> rdm;

      //! the actual orbital transform
      std::unique_ptr<simanneal::OrbitalTransform> orbtrans;

      //! the current unitary
      std::unique_ptr<simanneal::UnitaryMatrix> opt_unitary;

      //! our pseudo-random generator
      std::mt19937 mt;

      //! only rotation within these irreps (if not empty)
      std::vector<int> allow_irreps;
};

}

#endif /* LOCALMINIMIZER_H */

/*  vim: set ts=3 sw=3 expandtab :*/
