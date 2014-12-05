#ifndef OPTINDEX_H
#define OPTINDEX_H

#include <vector>

namespace CheMPS2 { class Hamiltonian; }

namespace simanneal
{

class OptIndex
{
    public:
        OptIndex(const int L, const int Group, const int * NORBin);
        OptIndex(const CheMPS2::Hamiltonian &);
        int getL() const; 
        int getNirreps() const;
        int getNORB(const int irrep) const;
        int getNstart(const int irrep) const;
        int * get_irrep_each_orbital();
        void Print() const;
    private:
        int Nirreps;
        int L;
        std::vector<int> NORB;
        std::vector<int> NORBcumulative;
        std::vector<int> irrep_each_orbital;
};

}

#endif /* OPTINDEX_H */
