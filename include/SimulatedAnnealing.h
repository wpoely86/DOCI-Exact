#ifndef SIM_ANNEAL_H
#define SIM_ANNEAL_H

#include <random>

#include "SymMolecule.h"
#include "DOCIHamtilonian.h"

namespace simanneal {
class OrbitalTransform;
class UnitaryMatrix;
}

namespace doci { class SimulatedAnnealing; }

class doci::SimulatedAnnealing
{
   public:
      SimulatedAnnealing(doci::Sym_Molecule &);

      SimulatedAnnealing(doci::Sym_Molecule &&);

      virtual ~SimulatedAnnealing();

      bool accept_function(double);

      void optimize();

      double calc_new_energy();

      void calc_energy();

      double get_energy() const;

      void Set_max_angle(double);

      void Set_delta_angle(double);

      void Set_start_temp(double);

      void Set_delta_temp(double);

      doci::DOCIHamiltonian& getHam() const;

      doci::Sym_Molecule& getMol() const;

      simanneal::OrbitalTransform& getOrbitaltf() const;

   private:

      //! Holds the current hamiltonian
      std::unique_ptr<doci::DOCIHamiltonian> ham;

      //! the actual orbital transform
      std::unique_ptr<simanneal::OrbitalTransform> orbtrans;

      //! the current unitary
      std::unique_ptr<simanneal::UnitaryMatrix> opt_unitary;

      //! energy of current iteration
      double energy;
      //! the start temperatur
      double start_temp;
      //! the change in temperatur between steps
      double delta_temp;
      //! the change in the angle allows in steps
      double delta_angle;
      //! the maximum allows angle
      double max_angle;
      //! the number of steps done
      unsigned int steps;
      //! max number of steps
      unsigned int max_steps;

      double cur_temp;

      //! the real random input (hopefully)
      std::random_device rd;
      //! our pseudo-random generator
      std::mt19937_64 mt;
};


#endif /* SIM_ANNEAL_H */

/* vim: set ts=3 sw=3 expandtab :*/
