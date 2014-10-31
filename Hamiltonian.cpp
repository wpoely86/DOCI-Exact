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
   mat.reset(new SparseMatrix_CRS(dim));
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
   mat.reset(new SparseMatrix_CRS(dim));
}

DOCIHamiltonian::DOCIHamiltonian(const DOCIHamiltonian &orig)
{
   permutations.reset(new Permutation(*orig.permutations));
   molecule.reset(new Molecule(*orig.molecule));
   mat.reset(new SparseMatrix_CRS(*orig.mat));
}

DOCIHamiltonian& DOCIHamiltonian::operator=(const DOCIHamiltonian &orig)
{
   permutations.reset(new Permutation(*orig.permutations));
   molecule.reset(new Molecule(*orig.molecule));
   mat.reset(new SparseMatrix_CRS(*orig.mat));

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
   return mat->gn();
}

/**
 * Build the (sparse) DOCIHamiltonian
 */
void DOCIHamiltonian::Build()
{
   auto &perm_bra = *permutations;
   perm_bra.reset();

   for(unsigned int i=0;i<mat->gn();++i)
   {
      const auto bra = perm_bra.get();

      mat->NewRow();

      // do all diagonal terms
      auto cur = bra;
      double tmp = 0;

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
         tmp += 2 * molecule->getT(s, s);

         // TEI: part a \bar a ; a \bar a
         tmp += molecule->getV(s, s, s, s);

         auto cur2 = cur; 

         while(cur2)
         {
            // select rightmost up state in the ket
            auto ksp2 = cur2 & (~cur2 + 1);
            // set it to zero
            cur2 ^= ksp2;

            // number of the orbital
            auto r = CountBits(ksp2-1);

            // s < r !! (avoid double counting)

            // TEI:
            // - a b ; a b
            // - a \bar b ; a \bar b
            // - \bar a \bar b ; \bar a \bar b
            // with a < b
            // The second term (ab|V|ba) is not possible in the second
            // case, so only a prefactor of 2 instead of 4.
            tmp += 4 * molecule->getV(r, s, r, s);
            tmp -= 2 * molecule->getV(r, s, s, r);
         }
      }

      mat->PushToRowNext(i, tmp);

      Permutation perm_ket(perm_bra);

      for(unsigned int j=i+1;j<mat->gn();++j)
      {
         const auto ket = perm_ket.next();

         const auto diff = bra ^ ket;

         // this means 4 orbitals are different
         if(CountBits(diff) == 2)
         {
            auto diff_c = diff;

            // select rightmost up state in the ket
            auto ksp1 = diff_c & (~diff_c + 1);
            // set it to zero
            diff_c ^= ksp1;

            auto ksp2 = diff_c & (~diff_c + 1);

            // number of the orbital
            auto r = CountBits(ksp1-1);
            auto s = CountBits(ksp2-1);

            // TEI: a \bar a ; b \bar b
            mat->PushToRowNext(j, molecule->getV(s, s, r, r));
         } 
      }

      perm_bra.next();
   }

   // closed the matrix properly
   mat->NewRow();
}


std::pair< double,std::vector<double> > DOCIHamiltonian::Diagonalize() const
{
   // dimension of the matrix
   int n = mat->gn();

   // number of eigenvalues to calculate
   int nev = 1;

   // reverse communication parameter, must be zero on first iteration
   int ido = 0;
   // standard eigenvalue problem A*x=lambda*x
   char bmat = 'I';
   // calculate the smallest algebraic eigenvalue
   char which[] = {'S','A'};
   // calculate until machine precision
   double tol = 0;

   // the residual vector
   std::unique_ptr<double []> resid(new double[n]);

   // the number of columns in v: the number of lanczos vector
   // generated at each iteration, ncv <= n
   // We use the answer to life, the universe and everything, if possible
   int ncv = 42;

   if( n < ncv )
      ncv = n;

   // v containts the lanczos basis vectors
   auto ldv = n;
   std::unique_ptr<double []> v(new double[ldv*ncv]);

   std::unique_ptr<int []> iparam(new int[11]);
   iparam[0] = 1;   // Specifies the shift strategy (1->exact)
   iparam[2] = 3*n; // Maximum number of iterations
   iparam[6] = 1;   /* Sets the mode of dsaupd.
                       1 is exact shifting,
                       2 is user-supplied shifts,
                       3 is shift-invert mode,
                       4 is buckling mode,
                       5 is Cayley mode. */

   std::unique_ptr<int []> ipntr(new int[11]); /* Indicates the locations in the work array workd
                                                  where the input and output vectors in the
                                                  callback routine are located. */

   // array used for reverse communication
   std::unique_ptr<double []> workd(new double[3*n]);
   for(int i=0;i<3*n;i++)
      workd[i] = 0;

   auto lworkl = ncv*(ncv+8); /* Length of the workl array */
   std::unique_ptr<double []> workl(new double[lworkl]);

   // info = 0: random start vector is used
   int info = 0; /* Passes convergence information out of the iteration
                    routine. */

   // rvec == 0 : calculate only eigenvalue
   // rvec > 0 : calculate eigenvalue and eigenvector
   int rvec = 0;

   // how many eigenvectors to calculate: 'A' => nev eigenvectors
   char howmny = 'A';

   std::unique_ptr<int []> select;
   // when howmny == 'A', this is used as workspace to reorder the eigenvectors
   if( howmny == 'A' )
      select.reset(new int[ncv]);

   // This vector will return the eigenvalues from the second routine, dseupd.
   std::unique_ptr<double []> d(new double[nev]);

   std::vector<double> eigv;

   if(rvec)
      eigv.resize(n);

   // not used if iparam[6] == 1
   double sigma;

   // first iteration
   dsaupd_(&ido, &bmat, &n, &which[0], &nev, &tol, resid.get(), &ncv, v.get(), &ldv, iparam.get(), ipntr.get(), workd.get(), workl.get(), &lworkl, &info);

   while( ido != 99 )
   {
      // matrix-vector multiplication
      mat->mvprod(workd.get()+ipntr[0]-1, workd.get()+ipntr[1]-1,0);

      dsaupd_(&ido, &bmat, &n, &which[0], &nev, &tol, resid.get(), &ncv, v.get(), &ldv, iparam.get(), ipntr.get(), workd.get(), workl.get(), &lworkl, &info);
   }

   if( info < 0 )
      std::cerr << "Error with dsaupd, info = " << info << std::endl;
   else if ( info == 1 )
      std::cerr << "Maximum number of Lanczos iterations reached." << std::endl;
   else if ( info == 3 )
      std::cerr << "No shifts could be applied during implicit Arnoldi update, try increasing NCV." << std::endl;

   dseupd_(&rvec, &howmny, select.get(), d.get(), eigv.data(), &ldv, &sigma, &bmat, &n, which, &nev, &tol, resid.get(), &ncv, v.get(), &ldv, iparam.get(), ipntr.get(), workd.get(), workl.get(), &lworkl, &info);

   if ( info != 0 )
      std::cerr << "Error with dseupd, info = " << info << std::endl;

   return std::make_pair(d[0], std::move(eigv));
}

/**
 * Convert the sparse matrix to a full matrix and use exact diagonalization
 * to find the eigenvalues and eigenvectors.
 * @return a pair of a vector with all the eigenvalues (sorted) and a matrix
 * with the corresponding orthonormal eigenvectors
 */
std::pair< std::vector<double>,helpers::matrix > DOCIHamiltonian::DiagonalizeFull() const
{
   char jobz = 'V';
   char uplo = 'U';
   int n = mat->gn();

   helpers::matrix fullmat(n, n);
   mat->ConvertToMatrix(fullmat);

   std::vector<double> eigs(n);

   int lwork = 3*n - 1;

   std::unique_ptr<double []> work (new double [lwork]);

   int info = 0;

   dsyev_(&jobz,&uplo,&n,fullmat.getpointer(),&n,eigs.data(),work.get(),&lwork,&info);

   if(info)
      std::cerr << "dsyev failed. info = " << info << std::endl;

   return std::make_pair(std::move(eigs), std::move(fullmat));
}

/**
 * Wrapper function for gcc buildin popcount
 * @param bits the unsigned long long to count the bits from
 * @return the number of ones in bits
 */
unsigned int DOCIHamiltonian::CountBits(mybitset bits)
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
int DOCIHamiltonian::CalcSign(unsigned int i,unsigned int j, const mybitset a)
{
    assert(i<j && "Order indices correctly!");

    // count the number of set bits between i and j in ket a
    auto sign = CountBits(( ((1<<j) - 1) ^ ((1<<(i+1)) - 1) ) & a);

    // when uneven, we get a minus sign
    if( sign & 0x1 )
        return -1;
    else
        return +1;
}

/* vim: set ts=3 sw=3 expandtab :*/
