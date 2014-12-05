#ifndef __OrbitalTransform__
#define __OrbitalTransform__

#include <assert.h>
#include <memory>
#include <vector>
#include "Irreps.h"
#include "OptIndex.h"

namespace CheMPS2 { class Hamiltonian; }

namespace simanneal {

class UnitaryMatrix;

class OrbitalTransform
{
    public :
        OrbitalTransform(const CheMPS2::Hamiltonian& ham);
        virtual ~OrbitalTransform() = default;

        void fillHamCI(CheMPS2::Hamiltonian& HamCI);	
        void fillConstAndTmat(CheMPS2::Hamiltonian& Ham) const;
        void buildOneBodyMatrixElements();
        void set_unitary(UnitaryMatrix& unit);
        CheMPS2::Hamiltonian& get_ham() { return (*_hamorig); }
        CheMPS2::Hamiltonian& get_ham() const { return (*_hamorig); }

        double TmatRotated(const int index1, const int index2) const;
        double get_norb(int irrep) const;
        double get_difference_orig(CheMPS2::Hamiltonian &hamin) const;
        void CheckDeviationFromUnitary() const;
        void update_unitary(double *step);

        UnitaryMatrix& get_unitary() { return (*_unitary); }
        UnitaryMatrix& get_unitary() const { return (*_unitary); }

    private:
        void rotate_old_to_new(std::unique_ptr<double []> * matrix);

        //! the orginal hamiltonian
        std::unique_ptr<CheMPS2::Hamiltonian> _hamorig;
        //! The rotation to perfrom on _hamorig to get the current hamiltonian
        std::unique_ptr<simanneal::UnitaryMatrix>  _unitary;

        OptIndex index;
        CheMPS2::Irreps SymmInfo;
        int numberOfIrreps;

        //! Some memory to do the one body work, allocated when needed
        std::vector< std::unique_ptr<double []> > QmatrixWork;
        // The one body matrix elements, used for rotations, only allocated when needed
        std::vector< std::unique_ptr<double []> > OneBodyMatrixElements;

        std::unique_ptr<double []> mem1;
        std::unique_ptr<double []> mem2;

};

}

#endif

/* vim: set ts=4 sw=4 expandtab :*/
