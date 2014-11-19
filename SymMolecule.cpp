#include <iostream>
#include <assert.h>

#include "SymMolecule.h"
// From CIFLOW 
#include <CIFLow/Hamiltonian.h>

doci::Sym_Molecule::Sym_Molecule(std::string filename) : ham(nullptr)
{
   ham.reset(new Hamiltonian(filename));
}

doci::Sym_Molecule::Sym_Molecule(const Sym_Molecule &orig)
{
   ham.reset(new Hamiltonian(*orig.ham));
}

// cannot do this in header as Hamiltonian is an incomplete type:
// https://stackoverflow.com/questions/13414652/forward-declaration-with-unique-ptr
doci::Sym_Molecule::~Sym_Molecule() = default;

doci::Sym_Molecule* doci::Sym_Molecule::clone() const
{
   return new doci::Sym_Molecule(*this);
}

double doci::Sym_Molecule::getT(int a, int b) const
{
   return ham->getTmat(a,b);
}

double doci::Sym_Molecule::getV(int a, int b, int c, int d) const
{
   return ham->getVmat(a,b,c,d);
}

double doci::Sym_Molecule::HF_Energy() const
{
   return ham->get_hf_energy();
}

double doci::Sym_Molecule::get_nucl_rep() const
{
   return ham->getEconst();
}

unsigned int doci::Sym_Molecule::get_n_sp() const
{
   return ham->getL();
}

unsigned int doci::Sym_Molecule::get_n_electrons() const
{
   return ham->getnup() + ham->getndown();
}

/**
 * Access operator to the actual Hamiltonian object
 * @return the Hamiltonian object
 */
Hamiltonian& doci::Sym_Molecule::getHamObject() const
{
   return *ham;
}

/* vim: set ts=3 sw=3 expandtab :*/
