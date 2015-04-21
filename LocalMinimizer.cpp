#include <cassert>
#include <chrono>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <hdf5.h>
#include <signal.h>
#include <cstring>
#include <sstream>

#include "LocalMinimizer.h"
#include "Hamiltonian.h"
#include "OptIndex.h"

/**
 * @param mol the molecular data to use
 */
doci::LocalMinimizer::LocalMinimizer(const doci::Sym_Molecule &mol)
{
   orbtrans.reset(new simanneal::OrbitalTransform(mol.getHamObject()));

   method.reset(new doci::DOCIHamiltonian(mol));

   rdm.reset(new doci::DM2(mol));

   energy = 0;

   conv_crit = 1e-6;
   conv_steps = 25;

   std::random_device rd;
   mt = std::mt19937(rd());

   // expect to be comma seperated list of allowed irreps
   char *irreps_env = getenv("v2DM_DOCI_ALLOWED_IRREPS");
   if(irreps_env && strlen(irreps_env) > 0)
   {
      std::string irreps_string = irreps_env;
      const std::string delim = ",";
      CheMPS2::Irreps syminfo(mol.getHamObject().getNGroup());

      auto start = 0U;
      auto end = irreps_string.find(delim);
      try
      { 
         while (true)
         {
            auto elem = irreps_string.substr(start, end - start);
            if(elem.empty())
               break;

            int cur_irrep = std::stoi(elem);

            if(cur_irrep >= 0 && cur_irrep < syminfo.getNumberOfIrreps())
               allow_irreps.push_back(cur_irrep);

            start = end + delim.length();

            if(end >= std::string::npos)
               break;

            end = irreps_string.find(delim, start);
         }
      } catch (std::exception& e) {
         std::cout << "Invalid value in v2DM_DOCI_ALLOWED_IRREPS" << std::endl;
      }

      std::sort(allow_irreps.begin(), allow_irreps.end());
      std::cout << "Allowed irreps: ";
      for(auto &elem: allow_irreps)
         std::cout << elem << " ";
      std::cout << std::endl;
   }
}

doci::LocalMinimizer::LocalMinimizer(doci::Sym_Molecule &&mol)
{
   method.reset(new doci::DOCIHamiltonian(mol));

   orbtrans.reset(new simanneal::OrbitalTransform(mol.getHamObject()));

   rdm.reset(new doci::DM2(mol));

   energy = 0;

   conv_crit = 1e-6;
   conv_steps = 50;

   std::random_device rd;
   mt = std::mt19937(rd());
}

doci::LocalMinimizer::~LocalMinimizer() = default;

/**
 * @return the real energy (calculated + nuclear repulsion)
 */
double doci::LocalMinimizer::get_energy() const
{
   return energy + method->getMolecule().get_nucl_rep();
}

/**
 * Calculate the energy with the current
 * molecular data
 */
void doci::LocalMinimizer::calc_energy()
{
   auto start = std::chrono::high_resolution_clock::now();
   method->Build();
   auto end = std::chrono::high_resolution_clock::now();

   std::cout << "Building took: " << std::fixed << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;

   start = std::chrono::high_resolution_clock::now();
   auto eig = method->Diagonalize();
   energy = eig.first;
   end = std::chrono::high_resolution_clock::now();
   std::cout << "E = " << eig.first + method->getMolecule().get_nucl_rep() << std::endl;

   std::cout << "Diagonalization took: " << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;

   auto perm = method->getPermutation();
   start = std::chrono::high_resolution_clock::now();
   rdm->Build(perm, eig.second);
   end = std::chrono::high_resolution_clock::now();

   std::cout << "Building 2DM took: " << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;
}

/**
 * Calculate the energy with the current
 * molecular data
 */
double doci::LocalMinimizer::calc_new_energy()
{
   auto *mol = static_cast<Sym_Molecule *> (&method->getMolecule());
   assert(mol && "Shit, NULL pointer");

   orbtrans->fillHamCI(mol->getHamObject());

   auto start = std::chrono::high_resolution_clock::now();
   method->Build();
   auto end = std::chrono::high_resolution_clock::now();

   std::cout << "Building took: " << std::fixed << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;

   start = std::chrono::high_resolution_clock::now();
   auto eig = method->Diagonalize();
   end = std::chrono::high_resolution_clock::now();
   std::cout << "E = " << eig.first + mol->get_nucl_rep() << std::endl;

   std::cout << "Diagonalization took: " << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;

   auto perm = method->getPermutation();
   start = std::chrono::high_resolution_clock::now();
   rdm->Build(perm, eig.second);
   end = std::chrono::high_resolution_clock::now();

   std::cout << "Building 2DM took: " << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;

   return eig.first;
}

double doci::LocalMinimizer::calc_new_energy(const doci::Sym_Molecule &new_ham)
{
   auto *mol = static_cast<Sym_Molecule *> (&method->getMolecule());
   assert(mol && "Shit, NULL pointer");

   mol->getHamObject() = new_ham.getHamObject();

   auto start = std::chrono::high_resolution_clock::now();
   method->Build();
   auto end = std::chrono::high_resolution_clock::now();

   std::cout << "Building took: " << std::fixed << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;

   start = std::chrono::high_resolution_clock::now();
   auto eig = method->Diagonalize();
   end = std::chrono::high_resolution_clock::now();
   std::cout << "E = " << eig.first + mol->get_nucl_rep() << std::endl;

   std::cout << "Diagonalization took: " << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;

   auto perm = method->getPermutation();
   start = std::chrono::high_resolution_clock::now();
   rdm->Build(perm, eig.second);
   end = std::chrono::high_resolution_clock::now();

   std::cout << "Building 2DM took: " << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;

   return eig.first;
}

simanneal::UnitaryMatrix& doci::LocalMinimizer::get_Optimal_Unitary()
{
   return orbtrans->get_unitary();
}

doci::Sym_Molecule& doci::LocalMinimizer::getHam() const
{
   auto *mol = static_cast<Sym_Molecule *> (&method->getMolecule());
   assert(mol && "Shit, NULL pointer");

   return *mol;
}

simanneal::OrbitalTransform& doci::LocalMinimizer::getOrbitalTf() const
{
   return *orbtrans;
}

std::vector< std::tuple<int,int,double,double> > doci::LocalMinimizer::scan_orbitals()
{
   auto start = std::chrono::high_resolution_clock::now();

   auto *mol = static_cast<Sym_Molecule *> (&method->getMolecule());
   assert(mol && "Shit, NULL pointer");

   const auto& ham2 = mol->getHamObject();
   std::function<double(int,int)> getT = [&ham2] (int a, int b) -> double { return ham2.getTmat(a,b); };
   std::function<double(int,int,int,int)> getV = [&ham2]  (int a, int b, int c, int d) -> double { return ham2.getVmat(a,b,c,d); };

   std::vector< std::tuple<int,int,double,double> > pos_rotations;
   // worst case: c1 symmetry
   pos_rotations.reserve(ham2.getL()*(ham2.getL()-1)/2);

   for(int k_in=0;k_in<ham2.getL();k_in++)
      for(int l_in=k_in+1;l_in<ham2.getL();l_in++)
         if(ham2.getOrbitalIrrep(k_in) == ham2.getOrbitalIrrep(l_in))
         {
            if(!allow_irreps.empty() && std::find(allow_irreps.begin(), allow_irreps.end(), ham2.getOrbitalIrrep(k_in)) == allow_irreps.end() )
               continue;

            auto found = rdm->find_min_angle(k_in,l_in,0.3,getT,getV);

            if(!found.second)
               // we hit a maximum
               found = rdm->find_min_angle(k_in,l_in,0.01,getT,getV);

            if(!found.second)
               // we're still stuck in a maximum, skip this!
               continue;

            // skip angles larger than Pi/2
            if(fabs(found.first)>M_PI/2.0)
               continue;

            double new_en = rdm->calc_rotate(k_in,l_in,found.first,getT,getV);

            assert(found.second && "Shit, maximum!");

            pos_rotations.push_back(std::make_tuple(k_in,l_in,found.first,new_en));
         }

   auto end = std::chrono::high_resolution_clock::now();

   std::cout << "Orbital scanning took: " << std::fixed << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;

   assert(pos_rotations.size()>0);

   return pos_rotations;
}

/**
 * Do the local minimization
 * @param dist_choice if set to true, we use choose_orbitals to choose
 * which pair of orbitals to use (instead of the lowest one)
 */
void doci::LocalMinimizer::Minimize(bool dist_choice)
{
   int converged = 0;
   double new_energy;

   // first run
   energy = calc_new_energy();

   auto start = std::chrono::high_resolution_clock::now();

   std::pair<int,int> prev_pair(0,0);

   int iters = 1;

   auto *mol = static_cast<Sym_Molecule *> (&method->getMolecule());
   assert(mol && "Shit, NULL pointer");

   auto& ham2 = mol->getHamObject();

   while(converged<conv_steps)
   {
      auto list_rots = scan_orbitals();

      std::sort(list_rots.begin(), list_rots.end(),
            [](const std::tuple<int,int,double,double> & a, const std::tuple<int,int,double,double> & b) -> bool
            {
            return std::get<3>(a) < std::get<3>(b);
            });

      for(auto& elem: list_rots)
         std::cout << std::get<0>(elem) << "\t" << std::get<1>(elem) << "\t" << std::get<3>(elem)+ham2.getEconst() << "\t" << std::get<2>(elem) << std::endl;

      int idx = 0;
      std::pair<int,int> tmp;

      if(dist_choice)
      {
         idx = choose_orbitalpair(list_rots);

         tmp = std::make_pair(std::get<0>(list_rots[idx]), std::get<1>(list_rots[idx]));

         if(tmp==prev_pair)
            idx = choose_orbitalpair(list_rots);

         tmp = std::make_pair(std::get<0>(list_rots[idx]), std::get<1>(list_rots[idx]));

         if(tmp==prev_pair)
            idx = 0;
      }

      tmp = std::make_pair(std::get<0>(list_rots[idx]), std::get<1>(list_rots[idx]));

      // don't do the same pair twice in a row
      if(tmp==prev_pair)
         idx++;

      const auto& new_rot = list_rots[idx];
      prev_pair = std::make_pair(std::get<0>(new_rot), std::get<1>(new_rot));

      if(dist_choice)
         std::cout << iters << " (" << converged << ") Chosen: " << idx << std::endl;

      assert(ham2.getOrbitalIrrep(std::get<0>(new_rot)) == ham2.getOrbitalIrrep(std::get<1>(new_rot)));
      // do Jacobi rotation twice: once for the Hamiltonian data and once for the Unitary Matrix
      orbtrans->DoJacobiRotation(ham2, std::get<0>(new_rot), std::get<1>(new_rot), std::get<2>(new_rot));
      orbtrans->get_unitary().jacobi_rotation(ham2.getOrbitalIrrep(std::get<0>(new_rot)), std::get<0>(new_rot), std::get<1>(new_rot), std::get<2>(new_rot));

      new_energy = calc_new_energy();

      std::stringstream h5_name;

      if(iters%10==0)
      {
         h5_name << getenv("SAVE_H5_PATH") << "/unitary-" << iters << ".h5";
         orbtrans->get_unitary().saveU(h5_name.str());
      }

      if(iters%25==0)
      {
         h5_name.str("");
         h5_name << getenv("SAVE_H5_PATH") << "/ham-" << iters << ".h5";
         ham2.save2(h5_name.str());

         h5_name.str("");
         h5_name << getenv("SAVE_H5_PATH") << "/rdm-" << iters << ".h5";
         rdm->WriteToFile(h5_name.str());
      }

      if(fabs(energy-new_energy)<conv_crit)
         converged++;

      std::cout << iters << " (" << converged << ")\tRotation between " << std::get<0>(new_rot) << "  " << std::get<1>(new_rot) << " over " << std::get<2>(new_rot) << " E_rot = " << std::get<3>(new_rot)+ham2.getEconst() << "  E = " << new_energy+ham2.getEconst() << "\t" << fabs(energy-new_energy) << std::endl;

      energy = new_energy;

      iters++;

      if(iters>1000)
      {
         std::cout << "Done 1000 steps, quiting..." << std::endl;
         break;
      }
   }

   auto end = std::chrono::high_resolution_clock::now();

   std::cout << "Minimization took: " << std::fixed << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;

   std::stringstream h5_name;
   h5_name << getenv("SAVE_H5_PATH") << "/optimale-uni.h5";
   get_Optimal_Unitary().saveU(h5_name.str());
}

double doci::LocalMinimizer::get_conv_crit() const
{
   return conv_crit;
}

void doci::LocalMinimizer::set_conv_crit(double crit)
{
   conv_crit = crit;
}

void doci::LocalMinimizer::set_conv_steps(int steps)
{
   conv_steps = steps;
}

/**
 * Choose a pair of orbitals to rotate over, according to the distribution of their relative
 * energy change.
 * @param orbs the list returned by scan_orbitals()
 * @return the index of the pair of orbitals in orbs
 */
int doci::LocalMinimizer::choose_orbitalpair(std::vector<std::tuple<int,int,double,double>> &orbs)
{
   std::uniform_real_distribution<double> dist(0, 1);

   const double choice = dist(mt);

   double norm = 0;

   for(auto &orb_pair: orbs)
      norm += (energy - std::get<3>(orb_pair));

   double cum = 0;
   for(int i=0;i<orbs.size();i++)
   {
      cum += (energy - std::get<3>(orbs[i]))/norm;
      if(choice < cum)
         return i;
   }

   assert(0 && "Should never ever be reached!");
   return -1;
}

const doci::DM2& doci::LocalMinimizer::get_DM2() const
{
   return *rdm;
}

/* vim: set ts=3 sw=3 expandtab :*/
