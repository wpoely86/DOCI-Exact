#include <stdexcept>
#include <assert.h>

#include "lapack.h"
#include "Hamiltonian.h"


/**
 * Constructor
 * @param perm the Permutation to use
 * @param mol the Molecule to use
 */
DOCIHamiltonian::DOCIHamiltonian(const Permutation &perm, const Molecule &mol)
{
   permutations.reset(new Permutation(perm));
   molecule.reset(new Molecule(mol));

   if(molecule->get_n_electrons() % 2 != 0)
      throw("We need even number of electrons!");

   auto dim = Permutation::CalcCombinations(molecule->get_n_sp(), molecule->get_n_electrons()/2);
   fullmat.reset(new helpers::matrix(dim, dim));
}

/**
 * Constructor.
 * We generate a Permutation object ourself based on the
 * Molecule in mol.
 * @param mol the Molecule to use
 */
DOCIHamiltonian::DOCIHamiltonian(const Molecule &mol)
{
   molecule.reset(new Molecule(mol));

   if(molecule->get_n_electrons() % 2 != 0)
      throw("We need even number of electrons!");

   permutations.reset(new Permutation(molecule->get_n_electrons()/2));

   auto dim = Permutation::CalcCombinations(molecule->get_n_sp(), molecule->get_n_electrons()/2);
   fullmat.reset(new helpers::matrix(dim, dim));
}

DOCIHamiltonian::DOCIHamiltonian(const DOCIHamiltonian &orig)
{
   permutations.reset(new Permutation(*orig.permutations));
   molecule.reset(new Molecule(*orig.molecule));
   mat.reset(new SparseMatrix_CCS(*orig.mat));
   fullmat.reset(new helpers::matrix(*orig.fullmat));
}

DOCIHamiltonian& DOCIHamiltonian::operator=(const DOCIHamiltonian &orig)
{
   permutations.reset(new Permutation(*orig.permutations));
   molecule.reset(new Molecule(*orig.molecule));
   mat.reset(new SparseMatrix_CCS(*orig.mat));
   fullmat.reset(new helpers::matrix(*orig.fullmat));

   return *this;
}

/**
 * @return the Molecule object
 */
Molecule const & DOCIHamiltonian::getMolecule() const
{
   return *molecule;
}

/**
 * @return the Permutation object
 */
Permutation const & DOCIHamiltonian::getPermutation() const
{
   return *permutations;
}

/**
 * @return the dimension of the hamiltonian matrix
 */
unsigned int DOCIHamiltonian::getdim() const
{
   return fullmat->getn();
}

/**
 * Build the (sparse) DOCIHamiltonian
 */
void DOCIHamiltonian::Build()
{
   auto &perm_bra = *permutations;

   for(unsigned int i=0;i<fullmat->getn();++i)
   {
      const auto bra = perm_bra.get();

      Permutation perm_ket(*permutations);

      for(unsigned int j=i;j<fullmat->getn();++j)
      {
         const auto ket = perm_ket.get();

         const auto diff = bra ^ ket;

         assert(CountBits(diff)%2 == 0);

         (*fullmat)(i,j) = 0;

         // this means 4 orbitals are different
         if(CountBits(diff) == 2)
         {

         }

         (*fullmat)(j,i) = (*fullmat)(i,j);

         perm_ket.next();
      }

      // do all diagonal terms
      auto cur = bra;

      // find occupied orbitals
      while(cur)
      {
         // select rightmost up state in the ket
         auto ksp = cur & (~cur + 1);
         // set it to zero
         cur ^= ksp;

         // number of the orbital
         auto s = CountBits(ksp-1);

         // OEI part
         (*fullmat)(i,i) += 2 * molecule->getT(s, s);

         // TEI: part a \bar a ; b \bar b (when a==b)
         (*fullmat)(i,i) += molecule->getV(s, s, s, s);

         auto cur2 = cur; 

         while(cur2)
         {
            // select rightmost up state in the ket
            auto ksp2 = cur2 & (~cur2 + 1);
            // set it to zero
            cur2 ^= ksp2;

            // number of the orbital
            auto r = CountBits(ksp2-1);

            assert(r>s);

            // TEI: part a b ; a b (4x)
            (*fullmat)(i,i) += molecule->getV(s, r, s, r) - molecule->getV(s, r, r, s);

            // TEI: part a \bar a ; b \bar b (when a!=b)
            (*fullmat)(i,i) += molecule->getV(s, s, r, r);
         }
      }

      perm_bra.next();
   }
}

std::vector<double> DOCIHamiltonian::diag()
{
   char jobz = 'N';
   char uplo = 'U';
   int n = fullmat->getn();

   std::vector<double> eigs(n);

   int lwork = 3*n - 1;

   std::unique_ptr<double []> work (new double [lwork]);

   int info = 0;

   dsyev_(&jobz,&uplo,&n,fullmat->getpointer(),&n,eigs.data(),work.get(),&lwork,&info);

   if(info)
      std::cerr << "dsyev failed. info = " << info << std::endl;

   return eigs;
}

/**
 * Wrapper function for gcc buildin popcount
 * @param bits the unsigned long long to count the bits from
 * @return the number of ones in bits
 */
unsigned int DOCIHamiltonian::CountBits(mybitset bits) const
{
#if defined(USELONG)
   return __builtin_popcountl(bits);
#elif defined(USELONGLONG)
   return __builtin_popcountll(bits);
#endif
}

/**
 * Calculate the sign by counting the number of set bits in a between position i and j
 * @param i the first position
 * @param j the second position, i < j
 * @param a the set of bits to calcalute on
 * @return the number of bits set between i and j in a
 */
int DOCIHamiltonian::CalcSign(unsigned int i,unsigned int j, const mybitset a) const
{
    assert(i<j);

    // count the number of set bits between i and j in ket a
    auto sign = CountBits(( ((1<<j) - 1) ^ ((1<<(i+1)) - 1) ) & a);

    // when uneven, we get a minus sign
    if( sign & 0x1 )
        return -1;
    else
        return +1;
}

/* vim: set ts=3 sw=3 expandtab :*/
