/* Copyright (C) 2014  Ward Poelmans

This file is part of DOCI-Exact.

Doci-Exact is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Hubbard-GPU is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with DOCI-Exact.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <getopt.h>

#include "Permutation.h"
#include "Molecule.h"
#include "DOCIHamtilonian.h"
#include "DM2.h"

// comment out if you don't need it
#include "SymMolecule.h"
#include "SimulatedAnnealing.h"
#include "Hamiltonian.h"
#include "OrbitalTransform.h"
#include "UnitaryMatrix.h"
#include "LocalMinimizer.h"

/**
 * @mainpage
 * This is an exact DOCI solver by means of a lanczos solver. We build the
 * Hamiltonian and store it as a sparse matrix (in Column Compressed Format).
 * The Permutation object generates the next basis state on the fly.
 * The Molecule object holds the atomic/molecular integrals to use.
 *
 * @author Ward Poelmans <wpoely86@gmail.com>
 * @version   0.6
 * @date      2014-2015
 * @copyright GNU Public License v3
 */

int main(int argc, char **argv)
{
    using namespace doci;
    using std::cout;
    using std::endl;

    cout.precision(10);

    std::string integralsfile = "mo-integrals.h5";
    std::string h5name = "rdm.h5";
    std::string unitary;
    bool simanneal = false;
    bool jacobirots = false;

    struct option long_options[] =
    {
        {"integrals",  required_argument, 0, 'i'},
        {"output",  required_argument, 0, 'o'},
        {"unitary",  required_argument, 0, 'u'},
        {"simulated-annealing",  no_argument, 0, 's'},
        {"jacobi-rotations",  no_argument, 0, 'j'},
        {"help",  no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int i,j;

    while( (j = getopt_long (argc, argv, "hi:o:su:j", long_options, &i)) != -1)
        switch(j)
        {
            case 'h':
            case '?':
                cout << "Usage: " << argv[0] << " [OPTIONS]\n"
                    "\n"
                    "    -i, --integrals=integrals-file  Set the input integrals file\n"
                    "    -o, --output=h5-file            Set the output filename for the RDM\n"
                    "    -s, --simulated-annealing       Use stimulated Annealing to find lowest energy\n"
                    "    -j, --jacobi-rotations          Use Jacobi Rotations to find lowest energy\n"
                    "    -u, --unitary                   Use this unitary to calc energy\n"
                    "    -h, --help                      Display this help\n"
                    "\n";
                return 0;
                break;
            case 'i':
                integralsfile = optarg;
                break;
            case 'o':
                h5name = optarg;
                break;
            case 's':
                simanneal = true;
                break;
            case 'j':
                jacobirots = true;
                break;
            case 'u':
                unitary = optarg;
                break;
        }

    if(simanneal && jacobirots)
    {
        cout << "Cannot do both simulated annealing and jacobi rotations, you have to choose!" << endl;

        return 2;
    }

    // make sure we have a save path, even if it's not specify already
    // This will not overwrite an already set SAVE_H5_PATH
    setenv("SAVE_H5_PATH", "./", 0);

    cout << "Reading: " << integralsfile << endl;

    if(!simanneal && !jacobirots)
    {
        Sym_Molecule mol(integralsfile);
        auto& ham_ints = mol.getHamObject();

        if(!unitary.empty())
        {
            cout << "Reading unitary " << unitary << endl;

            simanneal::OrbitalTransform orbtrans(ham_ints);

            orbtrans.get_unitary().loadU(unitary);
            orbtrans.fillHamCI(ham_ints);
        }

        DOCIHamiltonian ham(mol);

        auto start = std::chrono::high_resolution_clock::now();
        ham.Build();
        auto end = std::chrono::high_resolution_clock::now();

        cout << "Building took: " << std::fixed << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << endl;

        start = std::chrono::high_resolution_clock::now();
        auto eig2 = ham.Diagonalize();
        end = std::chrono::high_resolution_clock::now();

        cout << "Diagonalization took: " << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << endl;

        cout << "E = " << eig2.first + mol.get_nucl_rep() << endl;

        DM2 rdm(mol);
        auto perm = ham.getPermutation();
        start = std::chrono::high_resolution_clock::now();
        rdm.Build(perm, eig2.second);
        end = std::chrono::high_resolution_clock::now();

        cout << "Building 2DM took: " << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << endl;

        DM2 rdm_ham(mol);

        rdm_ham.BuildHamiltonian(mol);

        cout << "DM2 Energy = " << rdm.Dot(rdm_ham) + mol.get_nucl_rep() << endl;
        cout << "DM2 Trace = " << rdm.Trace() << endl;

        std::string h5_name = getenv("SAVE_H5_PATH");
        h5_name += "/rdm.h5";

        cout << "Writing 2DM to " << h5_name << endl;

        rdm.WriteToFile(h5_name);

        return 0;
    }

    if(simanneal)
    {
        Sym_Molecule mol(integralsfile);

        SimulatedAnnealing opt(mol);

        if(!unitary.empty())
        {
            cout << "Reading unitary " << unitary << endl;
            opt.getOrbitaltf().get_unitary().loadU(unitary);
        }

        opt.Set_start_temp(0.1);
        opt.Set_delta_temp(0.99);
        opt.Set_max_angle(1.3);
        opt.Set_delta_angle(0.999);

        auto start = std::chrono::high_resolution_clock::now();

        opt.optimize();

        auto end = std::chrono::high_resolution_clock::now();

        cout << "E = " << opt.get_energy() << endl;

        cout << "Optimization took: " << std::fixed << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << endl;

        opt.calc_new_energy();

        auto &ham = opt.getHam();
        auto &opt_mol = ham.getMolecule();

        auto eigs = ham.Diagonalize();
        cout << "E = " << eigs.first + mol.get_nucl_rep() << endl;

        DM2 rdm(opt_mol);
        auto perm = ham.getPermutation();
        start = std::chrono::high_resolution_clock::now();
        rdm.Build(perm, eigs.second);
        end = std::chrono::high_resolution_clock::now();

        cout << "Building 2DM took: " << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << endl;

        DM2 rdm_ham(opt_mol);
        rdm_ham.BuildHamiltonian(opt_mol);

        cout << "DM2 Energy = " << rdm.Dot(rdm_ham) + mol.get_nucl_rep() << endl;
        cout << "DM2 Trace = " << rdm.Trace() << endl;

        rdm.WriteToFile(h5name);
    }

    if(jacobirots)
    {
        Sym_Molecule mol(integralsfile);

        LocalMinimizer opt(mol);

        if(!unitary.empty())
        {
            cout << "Reading unitary " << unitary << endl;
            opt.getOrbitalTf().get_unitary().loadU(unitary);
        }

        opt.Minimize();

        cout << "E = " << opt.get_energy() << endl;

        DM2 rdm_ham(opt.getHam());
        rdm_ham.BuildHamiltonian(opt.getHam());
        const DM2 &rdm = opt.get_DM2();

        cout << "DM2 Energy = " << rdm.Dot(rdm_ham) + mol.get_nucl_rep() << endl;
        cout << "DM2 Trace = " << rdm.Trace() << endl;

        std::string h5_name = getenv("SAVE_H5_PATH");
        h5_name += "/rdm.h5";

        cout << "Writing 2DM to " << h5_name << endl;

        rdm.WriteToFile(h5name);
    }

/*     PSI_C1_Molecule mol(integralsfile);
 * 
 *     cout << "RHF energy = " << mol.HF_Energy() + mol.get_nucl_rep() << endl;
 * 
 *     DOCIHamiltonian ham(mol);
 * 
 *     auto start = std::chrono::high_resolution_clock::now();
 *     if(getenv("READ_SPARSE_H5_FILE"))
 *     {
 *         std::stringstream h5_name;
 *         h5_name << getenv("READ_SPARSE_H5_FILE");
 *         ham.ReadFromFile(h5_name.str());
 *     }
 *     else
 *         ham.Build();
 *     auto end = std::chrono::high_resolution_clock::now();
 * 
 *     if(getenv("SAVE_SPARSE_H5_FILE"))
 *     {
 *         std::stringstream h5_name;
 *         h5_name << getenv("SAVE_SPARSE_H5_FILE");
 *         ham.SaveToFile(h5_name.str());
 *     }
 * 
 *     cout << "Building took: " << std::fixed << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << endl;
 * 
 *     //    auto eig = ham.DiagonalizeFull();
 *     //
 *     //    for(unsigned int i=0;i<eig.first.size();i++)
 *     //        cout << i << "\t" << eig.first[i] + mol.get_nucl_rep() << endl;
 * 
 *     start = std::chrono::high_resolution_clock::now();
 *     auto eig2 = ham.Diagonalize();
 *     end = std::chrono::high_resolution_clock::now();
 * 
 *     cout << "Diagonalization took: " << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << endl;
 * 
 *     cout << "E = " << eig2.first + mol.get_nucl_rep() << endl;
 * 
 *     DM2 rdm(mol);
 *     auto perm = ham.getPermutation();
 *     start = std::chrono::high_resolution_clock::now();
 *     rdm.Build(perm, eig2.second);
 *     end = std::chrono::high_resolution_clock::now();
 * 
 *     cout << "Building 2DM took: " << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << endl;
 * 
 *     DM2 rdm_ham(mol);
 * 
 *     rdm_ham.BuildHamiltonian(mol);
 * 
 *     cout << "DM2 Energy = " << rdm.Dot(rdm_ham) + mol.get_nucl_rep() << endl;
 *     cout << "DM2 Trace = " << rdm.Trace() << endl;
 * 
 *     rdm.WriteToFile(h5name);
 */

    return 0;
}

/* vim: set ts=8 sw=4 tw=0 expandtab :*/
