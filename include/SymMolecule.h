#ifndef SYM_MOLECULE_H
#define SYM_MOLECULE_H

#include <iostream>
#include "Molecule.h"

// not part of our namespace
class Hamiltonian;

namespace doci {

/**
 * Molecular orbital matrix elements, with
 * point group symmetry (per irrep)
 * Comes from CheMPS2 -> CIFLOW
 * To generate input for this class, have a look at
 * https://github.com/SebWouters/CheMPS2/blob/master/mointegrals/mointegrals.cc_PRINT
 */
class Sym_Molecule: public Molecule
{
    public:
        Sym_Molecule(std::string filename);

        Sym_Molecule(const Sym_Molecule &);

        virtual ~Sym_Molecule();

        Sym_Molecule* clone() const;

        double getT(int, int) const;

        double getV(int, int, int, int) const;

        double HF_Energy() const;

        double get_nucl_rep() const;

        unsigned int get_n_sp() const;

        unsigned int get_n_electrons() const;

        Hamiltonian& getHamObject() const;

    private:

        std::unique_ptr<Hamiltonian> ham;
};

}

#endif /* SYM)MOLECULE_H */

/* vim: set ts=3 sw=3 expandtab :*/
