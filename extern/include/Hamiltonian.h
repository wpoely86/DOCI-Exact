/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013, 2014 Sebastian Wouters

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along
   with this program; if not, write to the Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef HAMILTONIAN_CHEMPS2_H
#define HAMILTONIAN_CHEMPS2_H

#include <string>
#include <memory>

#include "Irreps.h"
#include "TwoIndex.h"
#include "FourIndex.h"
#include "Options.h"

using std::string;

namespace simanneal { class OptIndex; }

namespace CheMPS2{
/** Hamiltonian class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date February 8, 2013
    
    Container class for the Hamiltonian matrix elements.
    
    \section ham_info Specific Hamiltonian information
    
    Class containing all Hamiltonian information:\n
    - L: the number of orbitals
    - groupNumber (in SymmInfo): the number of the Abelian point group symmetry with real-valued character table (see Irreps.h)
    - orb2irrep: array with the irrep number for each orbital
    - Econst: nuclear repulsion energy; or any constant part of the energy not contained in the 1- or 2-particle matrix elements
    - Tmat: 1-particle matrix elements; Tmat\f$_{a,b}\f$ = 0 if \f$I_a\f$ is different from \f$I_b\f$
    - Vmat: 2-particle matrix elements; Vmat\f$_{a,b,c,d}\f$ = 0 if \f$I_a \otimes I_b\f$ is not equal to \f$I_c \otimes I_d\f$; the matrix elements are not antisymmetrized and are stored with the convention that both (a & c) and (b & d) have the same spatial variable for the nuclear repulsion integral (physics notation).
    
    The targeted spin, particle number and point group symmetry are not defined here. For convenience, the second quantized formulation of the Hamiltonian is given here: \n
    \f$ \hat{H} = E_{const} + \sum\limits_{ij\sigma} T_{ij} \delta_{I_i,I_j} \hat{a}_{i \sigma}^{\dagger} \hat{a}_{j \sigma} + \frac{1}{2} \sum\limits_{ijkl\sigma\tau} V_{ijkl} \delta_{I_i \otimes I_j \otimes I_k \otimes I_l, I_{trivial}} \hat{a}_{i \sigma}^{\dagger} \hat{a}_{j \tau}^{\dagger} \hat{a}_{l \tau} \hat{a}_{k \sigma} \f$\n
    where the latin letters denote site-indices and the greek letters spin projections. This Hamiltonian preserves spin, spin projection, particle number, and Abelian point group symmetry (if its character table is real at least).
 */
   class Hamiltonian{

      friend class simanneal::OptIndex;

      public:
      
         //! Constructor
         /** \param Norbitals The number of orbitals (L)
             \param nGroup The group number
             \param OrbIrreps Pointer to array containing the orbital irreps */
         Hamiltonian(const int Norbitals, const int nGroup, const int * OrbIrreps);
         
         //! Constructor which loads a Psi4 text dump of the Hamiltonian from disk. A Psi4 text dump can be generated with the plugin psi4plugins/mointegrals.cc_PRINT.
         /** \param file_psi4text The filename of the Psi4 text dump of the Hamiltonian */
         Hamiltonian(const string file_psi4text);
         
         //! Constructor which loads a Hamiltonian from disk which was previously saved as a Psi4 text dump or in HDF5 format. A Psi4 text dump can be generated with the plugin psi4plugins/mointegrals.cc_PRINT. An HDF5 dump can be generated with the plugin psi4plugins/mointegrals.cc_SAVEHAM; or by (1) creating a Hamiltonian with one of the other constructors, (2) filling it with setEconst(), setTmat() and setVmat(), and (3) calling save().
         /** \param fileh5 If true, attempt to load a Hamiltonian in HDF5 format. All three filenames should be set then! If false, attempt to load a Hamiltonian which was previously saved as a Psi4 text dump. Only the first filename should be set then!
             \param main_file If fileh5, the HDF5 Hamiltonian parent filename. If not fileh5, the filename of the Psi4 text dump of the Hamiltonian.
             \param file_tmat The HDF5 Hamiltonian Tmat filename
             \param file_vmat The HDF5 Hamiltonian Vmat filename */
         Hamiltonian(const bool fileh5, const string main_file=HAMILTONIAN_ParentStorageName, const string file_tmat=HAMILTONIAN_TmatStorageName, const string file_vmat=HAMILTONIAN_VmatStorageName);

         Hamiltonian(const Hamiltonian &);

         Hamiltonian(Hamiltonian &&) = default;

         Hamiltonian& operator=(const Hamiltonian &);

         Hamiltonian& operator=(Hamiltonian &&) = default;

         //! Destructor
         virtual ~Hamiltonian() = default;
         
         //! Get the number of orbitals
         /** \return The number of orbitals */
         int getL() const;

         //! Get the group number
         /** \return The group number */
         int getNGroup() const;

         //! Get an orbital irrep number
         /** \param nOrb The orbital number
             \return The irrep of orbital nOrb */
         int getOrbitalIrrep(const int nOrb) const;
         
         //! Set the constant energy
         /** \param val The new constant energy */
         void setEconst(const double val);
         
         //! Set a Tmat element
         /** \param index1 The first index
             \param index2 The second index
             \param val The new Tmat element */
         void setTmat(const int index1, const int index2, const double val);
         
         //! Set a Vmat element
         /** \param index1 The first index
             \param index2 The second index
             \param index3 The third index
             \param index4 The fourth index
             \param val The new Vmat element */
         void setVmat(const int index1, const int index2, const int index3, const int index4, const double val);
         
         //! Add to Vmat element
         /** \param index1 The first index
             \param index2 The second index
             \param index3 The third index
             \param index4 The fourth index
             \param val The value which should be added */
         void addToVmat(const int index1, const int index2, const int index3, const int index4, const double val);
         
         //! Get the constant energy
         /** \return The constant part of the Hamiltonian (nuclear repulsion & condensed orbitals) */
         double getEconst() const;
         
         //! Get a Tmat element
         /** \param index1 The first index
             \param index2 The second index
             \return \f$T_{index1,index2}\f$ */
         double getTmat(const int index1, const int index2) const;
         
         //! Get a Vmat element
         /** \param index1 The first index
             \param index2 The second index
             \param index3 The third index
             \param index4 The fourth index
             \return \f$V_{index1,index2,index3,index4}\f$ */
         double getVmat(const int index1, const int index2, const int index3, const int index4) const;
         
         //! Save the Hamiltonian
         /** \param file_parent The HDF5 Hamiltonian parent filename
             \param file_tmat The HDF5 Hamiltonian Tmat filename
             \param file_vmat The HDF5 Hamiltonian Vmat filename */
         void save(const string file_parent=HAMILTONIAN_ParentStorageName, const string file_tmat=HAMILTONIAN_TmatStorageName, const string file_vmat=HAMILTONIAN_VmatStorageName) const;
         
         //! Load the Hamiltonian
         /** \param file_parent The HDF5 Hamiltonian parent filename
             \param file_tmat The HDF5 Hamiltonian Tmat filename
             \param file_vmat The HDF5 Hamiltonian Vmat filename */
         void read(const string file_parent=HAMILTONIAN_ParentStorageName, const string file_tmat=HAMILTONIAN_TmatStorageName, const string file_vmat=HAMILTONIAN_VmatStorageName);

         void save2(const string filename) const;

         void read2(const string filename);

         void setNe(int);

         int getNe() const;
         
         //! Debug check certain elements and sums
         void debugcheck() const;

         static Hamiltonian CreateFromH5(const string filename);

         //! set everything to zero
         void reset();
      
      private:
      
         //number of orbitals
         int L;

         //number of electrons
         int Ne;
         
         //symmetry info
         Irreps SymmInfo;
         
         //irrep of each orbital
         std::unique_ptr<int []> orb2irrep;
         
         //number of orbitals per irrep
         std::unique_ptr<int []> irrep2num_orb;
         
         //index of an orbital within irrep block
         std::unique_ptr<int []> orb2indexSy;
         
         //1-particle matrix elements
         std::unique_ptr<TwoIndex> Tmat;
         
         //2-particle matrix elements
         std::unique_ptr<FourIndex> Vmat;
         
         //Constant part of the Hamiltonian
         double Econst;
         
         //If filename=="LOADH5" in Hamiltonian::Hamiltonian then the HDF5 Hamiltonian is loaded
         void CreateAndFillFromH5(const string file_parent, const string file_tmat, const string file_vmat);
         
         //If filename!="LOADH5" in Hamiltonian::Hamiltonian then a Psi4 dump in the file with name "filename" is loaded
         void CreateAndFillFromPsi4dump(const string filename);
         
   };
}

#endif
