#include <iostream>
#include <assert.h>

#include "OptIndex.h"
#include "Hamiltonian.h"

using simanneal::OptIndex;

OptIndex::OptIndex(const int L, const int _Group, const int * NORBin)
{
    this->L = L;
    CheMPS2::Irreps SymmInfo;
    SymmInfo.setGroup(_Group);
    this->Nirreps = SymmInfo.getNumberOfIrreps();

    NORB.resize(Nirreps);
    NORBcumulative.resize(Nirreps+1);

    int sum_check = 0;
    NORBcumulative[0] = 0;
    for (int irrep=0; irrep<Nirreps; irrep++)
    {
        NORB[irrep] = NORBin[irrep];

        sum_check += NORB[irrep];

        NORBcumulative[irrep+1] = NORBcumulative[irrep] + NORB[irrep];
    }

    if (sum_check != L)
        std::cerr << "OptIndex::OptIndex : Sum over all orbitals is not L." << std::endl;

    irrep_each_orbital.resize(NORBcumulative[Nirreps]);

    for (int irrep=0; irrep<Nirreps; irrep++)
        for (int cnt=0; cnt<NORB[irrep]; cnt++)
            irrep_each_orbital[ NORBcumulative[irrep] + cnt ] = irrep;
}

OptIndex::OptIndex(const CheMPS2::Hamiltonian &ham)
{
    this->L = ham.L;
    this->Nirreps = ham.SymmInfo.getNumberOfIrreps();

    NORB.resize(Nirreps);
    NORBcumulative.resize(Nirreps+1);

    int sum_check = 0;
    NORBcumulative[0] = 0;
    for (int irrep=0; irrep<Nirreps; irrep++)
    {
        NORB[irrep] = ham.irrep2num_orb[irrep];

        sum_check += NORB[irrep];

        NORBcumulative[irrep+1] = NORBcumulative[irrep] + NORB[irrep];
    }

    if (sum_check != L)
        std::cerr << "OptIndex::OptIndex : Sum over all orbitals is not L." << std::endl;

    irrep_each_orbital.resize(NORBcumulative[Nirreps]);

    for (int irrep=0; irrep<Nirreps; irrep++)
        for (int cnt=0; cnt<NORB[irrep]; cnt++)
            irrep_each_orbital[ NORBcumulative[irrep] + cnt ] = irrep;
}

int OptIndex::getL() const { return L; }

int OptIndex::getNirreps() const { return Nirreps; }

int OptIndex::getNORB(const int irrep) const { return NORB[irrep]; }

int OptIndex::getNstart(const int irrep) const { return NORBcumulative[irrep]; }

int * OptIndex::get_irrep_each_orbital() { return irrep_each_orbital.data(); }

void OptIndex::Print() const
{
    using std::cout;
    using std::endl;

    cout << "NORB  = [ ";
    for (int irrep=0; irrep<Nirreps-1; irrep++)
        cout << NORB[irrep] << " , ";
    cout << NORB[Nirreps-1] << " ]" << endl;
}
