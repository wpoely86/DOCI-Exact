#include <stdexcept>
#include <sstream>
#include <chrono>
#include <omp.h>
#include <assert.h>

#include "lapack.h"
#include "DOCIHamtilonian.h"

using namespace doci;

/**
 * Constructor
 * @param perm the Permutation to use
 * @param mol the Molecule to use
 */
DOCIHamiltonian::DOCIHamiltonian(const Permutation &perm, const Molecule &mol)
{
   permutations.reset(new Permutation(perm));
   molecule.reset(mol.clone());

   if(molecule->get_n_electrons() % 2 != 0)
      throw("We need even number of electrons!");

   auto dim = Permutation::CalcCombinations(molecule->get_n_sp(), molecule->get_n_electrons()/2);
   mat.reset(new helpers::SparseMatrix_CRS(dim));
}

/**
 * Constructor.
 * We generate a Permutation object ourself based on the
 * Molecule in mol.
 * @param mol the Molecule to use
 */
DOCIHamiltonian::DOCIHamiltonian(const Molecule &mol)
{
   molecule.reset(mol.clone());

   if(molecule->get_n_electrons() % 2 != 0)
      throw("We need even number of electrons!");

   permutations.reset(new Permutation(molecule->get_n_electrons()/2));

   auto dim = Permutation::CalcCombinations(molecule->get_n_sp(), molecule->get_n_electrons()/2);
   mat.reset(new helpers::SparseMatrix_CRS(dim));
}

DOCIHamiltonian::DOCIHamiltonian(Molecule &&mol)
{
   molecule.reset(mol.move());

   if(molecule->get_n_electrons() % 2 != 0)
      throw("We need even number of electrons!");

   permutations.reset(new Permutation(molecule->get_n_electrons()/2));

   auto dim = Permutation::CalcCombinations(molecule->get_n_sp(), molecule->get_n_electrons()/2);
   mat.reset(new helpers::SparseMatrix_CRS(dim));
}


DOCIHamiltonian::DOCIHamiltonian(const DOCIHamiltonian &orig)
{
   permutations.reset(new Permutation(*orig.permutations));
   molecule.reset(orig.molecule->clone());
   mat.reset(new helpers::SparseMatrix_CRS(*orig.mat));
}

DOCIHamiltonian& DOCIHamiltonian::operator=(const DOCIHamiltonian &orig)
{
   permutations.reset(new Permutation(*orig.permutations));
   molecule.reset(orig.molecule->clone());
   mat.reset(new helpers::SparseMatrix_CRS(*orig.mat));

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
 * @return the Molecule object
 */
Molecule & DOCIHamiltonian::getMolecule()
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
   auto num_t = omp_get_max_threads();
   unsigned long long num_elems = (getdim()*1ul*(getdim()+1ul))/2;
   unsigned long long size_part = num_elems/num_t + 1;

   // every thread should process the lines between i and i+1 
   // with i the thread number
   std::vector<unsigned long long> workload(num_t+1);
   workload.front() = 0;
   workload.back() = getdim();

   for(int i=1;i<num_t;i++)
   {
      auto num_lines = workload[i-1];
      unsigned long long num_elems = 0;

      while(num_elems < size_part)
         num_elems += getdim() - num_lines++;

      if(num_lines > getdim())
         num_lines = getdim();

      workload[i] = num_lines;
   }

   std::vector< std::unique_ptr<helpers::SparseMatrix_CRS> > smat_parts(num_t);

   std::cout << "Running with " << num_t << " threads." << std::endl;

   permutations->reset();

#pragma omp parallel
   {
      auto start = std::chrono::high_resolution_clock::now();
      auto me = omp_get_thread_num();

      smat_parts[me].reset(new helpers::SparseMatrix_CRS(workload[me+1] - workload[me]));

      Permutation my_perm(*permutations);
      for(auto idx_begin=0;idx_begin<workload[me];++idx_begin)
         my_perm.next();

      // this costs some memory, make it an alias to
      // decrease memory usage but increase runtime
      auto my_mol = std::unique_ptr<Molecule> (molecule->clone());

      Build_iter(my_perm, (*smat_parts[me]), workload[me], workload[me+1], *my_mol);

      auto end = std::chrono::high_resolution_clock::now();

#pragma omp critical
      std::cout << me << "\t" << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;
   }

   mat->AddList(smat_parts);
}

/**
 * Internal method: this will iterate and build a part of the full sparse hamiltonian matrix
 * @param perm the start permutation to use
 * @param mat where to store the sparse matrix data
 * @param i_start the start point to iter
 * @param i_end the end point of the iterations
 * @param mol the molecule data to use
 */
void DOCIHamiltonian::Build_iter(Permutation &perm, helpers::SparseMatrix_CRS &mat,unsigned long long i_start, unsigned long long i_end, Molecule &mol)
{
   auto &perm_bra = perm;

   for(auto i=i_start;i<i_end;++i)
   {
      const auto bra = perm_bra.get();

      mat.NewRow();

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
         tmp += 2 * mol.getT(s, s);

         // TEI: part a \bar a ; a \bar a
         tmp += mol.getV(s, s, s, s);

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
            tmp += 4 * mol.getV(r, s, r, s);
            tmp -= 2 * mol.getV(r, s, s, r);
         }
      }

      mat.PushToRowNext(i, tmp);

      Permutation perm_ket(perm_bra);

      for(auto j=i+1;j<getdim();++j)
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
            mat.PushToRowNext(j, mol.getV(s, s, r, r));
         } 
      }

      perm_bra.next();
   }
}

/**
 * Calcalate the lowest eigenvalue and eigenvector using lanczos method.
 * We use arpack for this.
 * @return a pair of the lowest eigenvalue and corresponding normalized eigenvector
 */
std::pair< double,std::vector<double> > DOCIHamiltonian::Diagonalize() const
{
   double energy;
   std::vector<double> eigv(mat->gn());

   Diagonalize_arpack(energy,eigv,true);

   return std::make_pair(energy, std::move(eigv));
}

/**
 * Calcalate the lowest eigenvalue using lanczos method.
 * We use arpack for this.
 * @return the lowest eigenvalue
 */
double DOCIHamiltonian::CalcEnergy() const
{
   double energy;
   std::vector<double> eigv(0);

   Diagonalize_arpack(energy,eigv,false);

   return energy;
}

/**
 * Calcalate the lowest eigenvalue and (depending on eigvec) eigenvector using lanczos method.
 * We use arpack for this.
 * @param energy on return will hold the lowest eigenvalue
 * @param eigv on return will hold the lowest eigenvector
 * @param eigvec if true, calc the eigenvector and store in eigv
 */
void DOCIHamiltonian::Diagonalize_arpack(double &energy, std::vector<double> &eigv, bool eigvec) const
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
   if(eigvec)
      rvec = 1;

   // how many eigenvectors to calculate: 'A' => nev eigenvectors
   char howmny = 'A';

   std::unique_ptr<int []> select;
   // when howmny == 'A', this is used as workspace to reorder the eigenvectors
   if( howmny == 'A' )
      select.reset(new int[ncv]);

   // This vector will return the eigenvalues from the second routine, dseupd.
   std::unique_ptr<double []> d(new double[nev]);

   if(eigvec)
      eigv.resize(n);

   // not used if iparam[6] == 1
   double sigma;

   // first iteration
   dsaupd_(&ido, &bmat, &n, &which[0], &nev, &tol, resid.get(), &ncv, v.get(), &ldv, iparam.get(), ipntr.get(), workd.get(), workl.get(), &lworkl, &info);

   while( ido != 99 )
   {
      // matrix-vector multiplication
      mat->mvprod(workd.get()+ipntr[0]-1, workd.get()+ipntr[1]-1);

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

   energy = d[0];
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

   std::unique_ptr<helpers::matrix> fullmat(new helpers::matrix(n, n));
   mat->ConvertToMatrix(*fullmat);

   std::vector<double> eigs(n);

   int lwork = 3*n - 1;

   std::unique_ptr<double []> work (new double [lwork]);

   int info = 0;

   dsyev_(&jobz,&uplo,&n,fullmat->getpointer(),&n,eigs.data(),work.get(),&lwork,&info);

   if(info)
      std::cerr << "dsyev failed. info = " << info << std::endl;

   return std::make_pair(std::move(eigs), std::move(*fullmat));
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

/**
 * Save the sparse matrix to a file. Stores only the sparse matrix,
 * not number of electrons, ...
 * @param filename the name of the file
 */
void DOCIHamiltonian::SaveToFile(std::string filename) const
{
   mat->WriteToFile(filename.c_str(), "ham");
}

/**
 * Read the sparse matrix from a file. Reads only the sparse matrix,
 * not number of electrons, ...
 * You still need to construct a object first with a suited Molecule object.
 * @param filename the name of the file
 */
void DOCIHamiltonian::ReadFromFile(std::string filename)
{
   mat->ReadFromFile(filename.c_str(), "ham");
}

/**
 * Calculate the number lowest energy levels
 * @param number the number of energy levels to calculate
 * @return list of the energies
 */
std::vector<double> DOCIHamiltonian::CalcEnergy(int number) const
{
   // dimension of the matrix
   int n = mat->gn();

   // number of eigenvalues to calculate
   int nev = number;

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
   std::vector<double> d(nev);

   // not used if iparam[6] == 1
   double sigma;

   // first iteration
   dsaupd_(&ido, &bmat, &n, &which[0], &nev, &tol, resid.get(), &ncv, v.get(), &ldv, iparam.get(), ipntr.get(), workd.get(), workl.get(), &lworkl, &info);

   while( ido != 99 )
   {
      // matrix-vector multiplication
      mat->mvprod(workd.get()+ipntr[0]-1, workd.get()+ipntr[1]-1);

      dsaupd_(&ido, &bmat, &n, &which[0], &nev, &tol, resid.get(), &ncv, v.get(), &ldv, iparam.get(), ipntr.get(), workd.get(), workl.get(), &lworkl, &info);
   }

   if( info < 0 )
      std::cerr << "Error with dsaupd, info = " << info << std::endl;
   else if ( info == 1 )
      std::cerr << "Maximum number of Lanczos iterations reached." << std::endl;
   else if ( info == 3 )
      std::cerr << "No shifts could be applied during implicit Arnoldi update, try increasing NCV." << std::endl;

   dseupd_(&rvec, &howmny, select.get(), d.data(), 0, &ldv, &sigma, &bmat, &n, which, &nev, &tol, resid.get(), &ncv, v.get(), &ldv, iparam.get(), ipntr.get(), workd.get(), workl.get(), &lworkl, &info);

   if ( info != 0 )
      std::cerr << "Error with dseupd, info = " << info << std::endl;

   return d;
}

/* vim: set ts=3 sw=3 expandtab :*/
