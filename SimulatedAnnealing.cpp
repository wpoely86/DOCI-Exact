#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cassert>
#include <hdf5.h>

#include "SimulatedAnnealing.h"

#include "Hamiltonian.h"
#include "UnitaryMatrix.h"
#include "OrbitalTransform.h"

using namespace simanneal;

/**
 * You still need to set the max_angle, delta_angle, start_temp and
 * delta_temp after creating the object.
 * @param mol the molecular data to use
 */
doci::SimulatedAnnealing::SimulatedAnnealing(doci::Sym_Molecule &mol)
{
   ham.reset(new DOCIHamiltonian(mol));

   OptIndex index(mol.getHamObject());

   opt_unitary.reset(new UnitaryMatrix(index));

   orbtrans.reset(new OrbitalTransform(mol.getHamObject()));

   mt = std::mt19937_64(rd());

   steps = 0;
   energy = 0;
   max_steps = 20000;
}

doci::SimulatedAnnealing::SimulatedAnnealing(doci::Sym_Molecule &&mol)
{
   ham.reset(new DOCIHamiltonian(mol));

   OptIndex index(mol.getHamObject());

   opt_unitary.reset(new UnitaryMatrix(index));

   orbtrans.reset(new OrbitalTransform(mol.getHamObject()));

   mt = std::mt19937_64(rd());

   steps = 0;
   energy = 0;
   max_steps = 20000;
}

doci::SimulatedAnnealing::~SimulatedAnnealing() = default;

/**
 * @return the real energy (calculated + nuclear repulsion)
 */
double doci::SimulatedAnnealing::get_energy() const
{
   return energy + ham->getMolecule().get_nucl_rep();
}

void doci::SimulatedAnnealing::Set_max_angle(double max_angle)
{
   this->max_angle = max_angle;
}

void doci::SimulatedAnnealing::Set_delta_angle(double delta_angle)
{
   this->delta_angle = delta_angle;
}

void doci::SimulatedAnnealing::Set_start_temp(double start_temp)
{
   this->start_temp = start_temp;
}

void doci::SimulatedAnnealing::Set_delta_temp(double delta_temp)
{
   this->delta_temp = delta_temp;
}

/**
 * Decide wether or not to accept the new energy
 * @param e_new the new energy
 * @return accept or not
 */
bool doci::SimulatedAnnealing::accept_function(double e_new)
{
   std::uniform_real_distribution<double> dist_accept(0, 1);

   if(e_new < energy)
      return true;
   else
   {
      double chance = std::exp((energy - e_new) / cur_temp);

      if ( dist_accept(mt) * (1+chance) > chance)
         //biggest chance for failing
         return false;
      else 
         return true;
   }
}

/**
 * Calculate the energy with the current
 * molecular data
 */
void doci::SimulatedAnnealing::calc_energy()
{
   ham->Build();

   energy = ham->CalcEnergy();
}

/**
 * Calculate the energy with the current
 * molecular data
 */
double doci::SimulatedAnnealing::calc_new_energy()
{
   // cast from Molecule to Sym_Molecule
   auto *mol = static_cast<Sym_Molecule *> (&ham->getMolecule());
   assert(mol && "Shit, NULL pointer");

   orbtrans->fillHamCI(mol->getHamObject());

   ham->Build();

   return ham->CalcEnergy();
}

/**
 * Do the simulated annealing
 */
void doci::SimulatedAnnealing::optimize()
{
   std::uniform_int_distribution<int> dist(0, ham->getMolecule().get_n_sp()-1);
   std::uniform_real_distribution<double> dist_angles(0, 1);

   auto *mol = static_cast<Sym_Molecule *> (&ham->getMolecule());
   auto &ham_data = mol->getHamObject();

   unsigned int unaccepted = 0;

   energy = calc_new_energy();
   std::cout << "Starting energy = " << get_energy() << std::endl;

   double lowest_energy = energy;
   cur_temp = start_temp;

   helpers::matrix sample_pairs(ham_data.getL(), ham_data.getL());
   sample_pairs = 0;

   for(unsigned int i=0;i<ham_data.getL();i++)
      for(unsigned int j=i+1;j<ham_data.getL();j++)
         if(ham_data.getOrbitalIrrep(i) != ham_data.getOrbitalIrrep(j))
            sample_pairs(i,j) = sample_pairs(j,i) = -1;

   auto start = std::chrono::high_resolution_clock::now();

   unsigned int i;
   for(i=0;i<max_steps;i++)
   {
      auto orb1 = dist(mt);
      auto orb2 = dist(mt);

      if(orb1 != orb2 && ham_data.getOrbitalIrrep(orb1) == ham_data.getOrbitalIrrep(orb2))
      {
         sample_pairs(orb1, orb2) += 1;
         sample_pairs(orb2, orb1) += 1;

         // between -1 and 1 but higher probablity to be close to zero (seems to work better)
         auto cur_angle = max_angle * (dist_angles(mt) - dist_angles(mt));

         std::cout << i << "\tT=" << cur_temp << "\tOrb1=" << orb1 << "\tOrb2=" << orb2 << "  Over " << cur_angle << std::endl;

         orbtrans->get_unitary().jacobi_rotation(ham_data.getOrbitalIrrep(orb1), orb1, orb2, cur_angle);

         auto new_energy = calc_new_energy();

         if(new_energy < lowest_energy)
            lowest_energy = new_energy;

         std::cout << "T=" << cur_temp << "\tNew energy = " << new_energy + mol->get_nucl_rep() << "\t Old energy = " << get_energy();

         if(accept_function(new_energy))
         {
            energy = new_energy;
            std::cout << "\t=> Accepted" << std::endl;
         }
         else
         {
            unaccepted++;
            std::cout << "\t=> Unaccepted, " << unaccepted << std::endl;
            orbtrans->get_unitary().jacobi_rotation(ham_data.getOrbitalIrrep(orb1), orb1, orb2, -1*cur_angle);
         }

         cur_temp *= delta_temp;
         max_angle *= delta_angle;

         if(unaccepted > 1500)
         {
            std::cout << "Too many unaccepted, stopping" << std::endl;
            break;
         }

      }

   }

   auto end = std::chrono::high_resolution_clock::now();

   std::cout << "Bottom was " << lowest_energy + mol->get_nucl_rep() << std::endl;
   std::cout << "Final energy = " << get_energy() << std::endl;
   std::cout << "Sim anneal runtime: " << std::fixed << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;

   std::stringstream h5_name;
   if(getenv("SAVE_H5_PATH"))
      h5_name << getenv("SAVE_H5_PATH") << "/unitary-final-" << i << ".h5";
   else
      h5_name << "unitary-final-" << i << ".h5";
   orbtrans->get_unitary().saveU(h5_name.str());

   for(unsigned int i=0;i<ham_data.getL();i++)
      for(unsigned int j=i+1;j<ham_data.getL();j++)
         if(sample_pairs(i,j) >= 0)
            std::cout << ham_data.getOrbitalIrrep(i) << "\t" << i << "\t" << j << "\t" << sample_pairs(i,j) << std::endl;
}

doci::DOCIHamiltonian& doci::SimulatedAnnealing::getHam() const
{
   return *ham;
}

simanneal::OrbitalTransform& doci::SimulatedAnnealing::getOrbitaltf() const
{
   return *orbtrans;
}

/* vim: set ts=3 sw=3 expandtab :*/
