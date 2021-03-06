#include <algorithm>
#include <functional>
#include <hdf5.h>
#include <assert.h>
#include <omp.h>
#include <chrono>
#include <cmath>

#include "DM2.h"
#include "DOCIHamtilonian.h"
#include "lapack.h"

using namespace doci;

std::unique_ptr<helpers::matrix> doci::DM2::sp2tp = nullptr;
std::unique_ptr<helpers::matrix> doci::DM2::tp2sp = nullptr;

// this helps to check the return codes of HDF5 calls
#define HDF5_STATUS_CHECK(status) if(status < 0) std::cerr << __FILE__ << ":" << __LINE__ << ": Problem with writing to file. Status code=" << status << std::endl;

/**
 * Create second order density matrix for n_sp single particle levels
 * (with spin degeneracy), meaning: n_sp doci levels.
 * @param n_sp the number of doci levels
 * @param n the number of particles
 */
DM2::DM2(unsigned int n_sp, unsigned int n)
{
   if(!sp2tp || !tp2sp)
      fill_lists(n_sp);

   block.reset(new helpers::matrix(n_sp, n_sp));
   diag.resize((n_sp*(n_sp-1))/2);

   this->N =  n;
}


/**
 * Create second order density matrix for a DOCI 
 * wavefunction based on the Molecule in mol
 * @param mol the molecule to use
 */
DM2::DM2(const Molecule &mol)
{
   auto n_sp = mol.get_n_sp();

   if(!sp2tp || !tp2sp)
      fill_lists(n_sp);

   block.reset(new helpers::matrix(n_sp, n_sp));
   diag.resize((n_sp*(n_sp-1))/2);
   N = mol.get_n_electrons();
}

DM2::DM2(const DM2 &orig)
{
   assert(sp2tp && tp2sp);

   block.reset(new helpers::matrix(*orig.block));
   diag = orig.diag;
   N = orig.N;
}

DM2::DM2(DM2 &&orig)
{
   block = std::move(orig.block);
   diag = std::move(orig.diag);
   N = orig.N;
}

DM2& DM2::operator=(const DM2 &orig)
{
   block.reset(new helpers::matrix(*orig.block));
   diag = orig.diag;
   N = orig.N;

   return *this;
}

DM2& DM2::operator=(DM2 &&orig)
{
   block = std::move(orig.block);
   diag = std::move(orig.diag);
   N = orig.N;

   return *this;
}

/**
 * Put all elements equal to the same value
 * @param val the value to use
 * @return *this
 */
DM2& DM2::operator=(double val)
{ 
   (*block) = val;
   std::fill(diag.begin(), diag.end(), val);

   return *this;
}

/**
 * Access the element \f$\hat a^+_a \hat a^+_b \hat a_d \hat a_c\f$
 * of the second order density matrix
 * @param a the first sp index
 * @param b the second sp index
 * @param c the thirth sp index
 * @param d the fourth sp index
 * @return the corresponding DM2 value
 */
double DM2::operator()(int a, int b, int c, int d) const
{
   if(a==b || c==d)
      return 0;

   int sign = 1;

   if(a > b)
      sign *= -1;

   if(c > d)
      sign *= -1;

   int i = (*sp2tp)(a,b);
   int j = (*sp2tp)(c,d);

   if(i<block->getn() && j<block->getn())
      return sign * (*block)(i,j);
   else if(i==j)
      return sign * diag[(i-block->getn())%diag.size()];
   else
      return 0;
}

/**
 * Add two DM2 object elementwise
 * @param a the other DM2 object
 * @return the sum of *this and a
 */
DM2& DM2::operator+=(const DM2 &a)
{
   (*block) += (*a.block);
//   std::transform(diag.begin(), diag.end(), a.diag.begin(), diag.begin(), std::plus<double>());

   int dim = diag.size();
   int inc = 1;
   double alpha = 1.0;
   double *tmp = const_cast<double *>(a.diag.data());

   daxpy_(&dim,&alpha,tmp,&inc,diag.data(),&inc);

   return *this;
}

/**
 * Build the sp <-> tp lists
 * @param n_sp the number of DOCI levels
 */
void DM2::fill_lists(unsigned int n_sp)
{
   auto L = n_sp;
   auto M = 2*n_sp;
   auto n_tp = M*(M-1)/2;

   sp2tp.reset(new helpers::matrix(M,M));
   (*sp2tp) = -1; // if you use something you shouldn't, this will case havoc

   tp2sp.reset(new helpers::matrix(n_tp,2));
   (*tp2sp) = -1; // if you use something you shouldn't, this will case havoc

   auto tel = 0;

   // a \bar a
   for(int a=0;a<L;a++)
      (*sp2tp)(a,a+L) = (*sp2tp)(a+L,a) = tel++;

   // a b
   for(int a=0;a<L;a++)
      for(int b=a+1;b<L;b++)
         (*sp2tp)(a,b) = (*sp2tp)(b,a) = tel++;

   // \bar a \bar b
   for(int a=L;a<M;a++)
      for(int b=a+1;b<M;b++)
         (*sp2tp)(a,b) = (*sp2tp)(b,a) = tel++;

   // a \bar b ; a \bar b
   for(int a=0;a<L;a++)
      for(int b=L+a+1;b<M;b++)
         if(a%L!=b%L)
            (*sp2tp)(a,b) = (*sp2tp)(b,a) = tel++;

   // \bar a b ; \bar a b
   for(int a=L;a<M;a++)
      for(int b=a%L+1;b<L;b++)
         if(a%L!=b%L)
            (*sp2tp)(a,b) = (*sp2tp)(b,a) = tel++;

   assert(tel == n_tp);

   for(int a=0;a<M;a++)
      for(int b=a+1;b<M;b++)
      {
         (*tp2sp)((*sp2tp)(a,b),0) = a;
         (*tp2sp)((*sp2tp)(a,b),1) = b;
      }
}

/**
 * Write the current DM2 to a HDF5 file
 * @param filename the name of the file
 */
void DM2::WriteToFile(std::string filename) const
{
   hid_t       file_id, group_id, dataset_id, attribute_id, dataspace_id;
   herr_t      status;

   file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   HDF5_STATUS_CHECK(file_id);

   group_id = H5Gcreate(file_id, "/RDM", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   HDF5_STATUS_CHECK(group_id);

   hsize_t dimblock = block->getn() * block->getn();

   dataspace_id = H5Screate_simple(1, &dimblock, NULL);

   dataset_id = H5Dcreate(group_id, "Block", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   HDF5_STATUS_CHECK(dataspace_id);

   double *data = const_cast<helpers::matrix &>(*block).getpointer();

   status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
   HDF5_STATUS_CHECK(status);

   status = H5Sclose(dataspace_id);
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);

   dimblock = diag.size();

   dataspace_id = H5Screate_simple(1, &dimblock, NULL);

   dataset_id = H5Dcreate(group_id, "Vector", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   HDF5_STATUS_CHECK(dataspace_id);

   status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, diag.data());
   HDF5_STATUS_CHECK(status);

   status = H5Sclose(dataspace_id);
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);


   dataspace_id = H5Screate(H5S_SCALAR);
   unsigned int L = block->getn();

   attribute_id = H5Acreate (group_id, "L", H5T_STD_U32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Awrite (attribute_id, H5T_NATIVE_UINT, &L );
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   attribute_id = H5Acreate (group_id, "N", H5T_STD_U32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Awrite (attribute_id, H5T_NATIVE_UINT, &N );
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   status = H5Sclose(dataspace_id);
   HDF5_STATUS_CHECK(status);

   status = H5Gclose(group_id);
   HDF5_STATUS_CHECK(status);

   status = H5Fclose(file_id);
   HDF5_STATUS_CHECK(status);
}

/**
 * Read a DM2 in from a HDF5 file.
 * @param filename the file to use
 * @return a new DM2 object with the data from the HDF5 file
 */
DM2 DM2::ReadFromFile(std::string filename)
{
   hid_t       file_id, dataset_id, group_id, attribute_id;
   herr_t      status;

   file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   HDF5_STATUS_CHECK(file_id);

   group_id = H5Gopen(file_id, "/RDM", H5P_DEFAULT);
   HDF5_STATUS_CHECK(group_id);

   attribute_id = H5Aopen(group_id, "L", H5P_DEFAULT);
   HDF5_STATUS_CHECK(attribute_id);

   unsigned int L;
   status = H5Aread(attribute_id, H5T_NATIVE_UINT, &L);
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   attribute_id = H5Aopen(group_id, "N", H5P_DEFAULT);
   HDF5_STATUS_CHECK(attribute_id);

   unsigned int N;
   status = H5Aread(attribute_id, H5T_NATIVE_UINT, &N);
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   DM2 dm2(L,N);

   dataset_id = H5Dopen(file_id, "/RDM/Block", H5P_DEFAULT);
   HDF5_STATUS_CHECK(dataset_id);

   status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dm2.block->getpointer());
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);

   dataset_id = H5Dopen(file_id, "/RDM/Vector", H5P_DEFAULT);
   HDF5_STATUS_CHECK(dataset_id);

   status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dm2.diag.data());
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);

   status = H5Gclose(group_id);
   HDF5_STATUS_CHECK(status);

   status = H5Fclose(file_id);
   HDF5_STATUS_CHECK(status);

   return dm2;
}

/**
 * @return number of particles
 */
unsigned int DM2::get_n_electrons() const
{
   return N;
}

/**
 * @return number of DOCI levels
 */
unsigned int DM2::get_n_sp() const
{
   return block->getn();
}

/**
 * Build the second order density matrix from a DOCI wavefunction
 * using the Permutation object and the ground state eigen vector
 * @param perm the Permutation object to use
 * @param eigv the eigenvector to build the DM2 from
 */
void DM2::Build(Permutation &perm, std::vector<double> &eigv)
{
   auto num_t = omp_get_max_threads();
   unsigned long long num_elems = (eigv.size()*1ul*(eigv.size()+1ul))/2;
   unsigned long long size_part = num_elems/num_t + 1;

   // every thread should process the lines between i and i+1 
   // with i the thread number
   std::vector<unsigned int> workload(num_t+1);
   workload.front() = 0;
   workload.back() = eigv.size();

   for(int i=1;i<num_t;i++)
   {
      auto num_lines = workload[i-1];
      unsigned long long num_elems = 0;

      while(num_elems < size_part)
         num_elems += eigv.size() - num_lines++;

      if(num_lines > eigv.size())
         num_lines = eigv.size();

      workload[i] = num_lines;
   }

   std::vector< std::unique_ptr<DM2> > dm2_parts(num_t);

   std::cout << "Running with " << num_t << " threads." << std::endl;

   perm.reset();

#pragma omp parallel
   {
      auto start = std::chrono::high_resolution_clock::now();
      auto me = omp_get_thread_num();

      dm2_parts[me].reset(new DM2(block->getn(),N));
      (*dm2_parts[me]) = 0;

      Permutation my_perm(perm);

      for(auto idx_begin=0;idx_begin<workload[me];++idx_begin)
         my_perm.next();

      auto vec_copy = eigv;

      build_iter(my_perm, vec_copy, workload[me], workload[me+1], (*dm2_parts[me]));

      auto end = std::chrono::high_resolution_clock::now();

#pragma omp critical
      std::cout << me << "\t" << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;
   }

   // add everything
   (*this) = 0;
   for(auto &cur_dm2: dm2_parts)
      (*this) += (*cur_dm2);
}

void DM2::build_iter(Permutation& perm, std::vector<double> &eigv, unsigned int i_start, unsigned int i_end, DM2 &cur_2dm)
{
   auto& perm_bra = perm;

   for(unsigned int i=i_start;i<i_end;++i)
   {
      const auto bra = perm_bra.get();

      Permutation perm_ket(perm_bra);

      auto cur = bra;

      // find occupied orbitals
      while(cur)
      {
         // select rightmost up state in the ket
         auto ksp = cur & (~cur + 1);
         // set it to zero
         cur ^= ksp;

         // number of the orbital
         auto s = DOCIHamiltonian::CountBits(ksp-1);

         // in this case: s == (*sp2tp)(s,s+L)
         (*cur_2dm.block)(s,s) += eigv[i] * eigv[i];

         auto cur2 = cur; 

         while(cur2)
         {
            // select rightmost up state in the ket
            auto ksp2 = cur2 & (~cur2 + 1);
            // set it to zero
            cur2 ^= ksp2;

            // number of the orbital
            auto r = DOCIHamiltonian::CountBits(ksp2-1);

            unsigned int idx = (*sp2tp)(r,s);
            // find correct relative index
            idx -= block->getn();
            idx %= diag.size();

            cur_2dm.diag[idx] += eigv[i] * eigv[i];
         }
      }


      for(unsigned int j=i+1;j<eigv.size();++j)
      {
         const auto ket = perm_ket.next();

         const auto diff = bra ^ ket;

         // this means 4 orbitals are different
         if(DOCIHamiltonian::CountBits(diff) == 2)
         {
            auto diff_c = diff;

            // select rightmost up state in the ket
            auto ksp1 = diff_c & (~diff_c + 1);
            // set it to zero
            diff_c ^= ksp1;

            auto ksp2 = diff_c & (~diff_c + 1);

            // number of the orbital
            auto r = DOCIHamiltonian::CountBits(ksp1-1);
            auto s = DOCIHamiltonian::CountBits(ksp2-1);

            // in this case:
            // r == (*sp2tp)(r,r+L)
            // s == (*sp2tp)(s,s+L)
            (*cur_2dm.block)(r,s) += eigv[i] * eigv[j];
            (*cur_2dm.block)(s,r) += eigv[i] * eigv[j];
         }
      }

      perm_bra.next();
   }
}

std::ostream &operator<<(std::ostream &output,doci::DM2 &dm2)
{
   auto L = dm2.get_n_sp();
   output << "Block: " << std::endl;
   for(int i=0;i<L;i++)
      for(int j=i;j<L;j++)
         output << i << "\t" << j << "\t|\t" << (*dm2.tp2sp)(i,0) << "  " <<  (*dm2.tp2sp)(i,1) << " ; " <<  (*dm2.tp2sp)(j,0) << "  " <<  (*dm2.tp2sp)(j,1) << "\t\t" << (*dm2.block)(i,j) << std::endl;

   output << std::endl;

   output << "Vector (4x): " << std::endl;
   for(int i=0;i<dm2.diag.size();i++)
      output << i << "\t|\t" << (*dm2.tp2sp)(L+i,0) << "  " << (*dm2.tp2sp)(L+i,1) << "\t\t" << dm2.diag[i] << std::endl;

   return output;
}

/**
 * Build the reducted hamiltonian for Molecule mol
 * @param mol the molecule to use
 */
void DM2::BuildHamiltonian(const Molecule &mol)
{
   auto L = block->getn();

   // make our life easier
   auto calc_elem = [this,&mol,L] (int i, int j) {
      int a = (*tp2sp)(i,0);
      int b = (*tp2sp)(i,1);
      int c = (*tp2sp)(j,0);
      int d = (*tp2sp)(j,1);

      int a_ = a % L;
      int b_ = b % L;
      int c_ = c % L;
      int d_ = d % L;

      double result = 0;

      // sp terms
      if(i==j)
         result += (mol.getT(a_,a_) + mol.getT(b_,b_))/(N - 1.0);

      // tp terms:

      // a \bar a ; b \bar b
      if(b==(a+L) && d==(c+L))
         result += mol.getV(a_,b_,c_,d_);

      // a b ; a b
      // \bar a \bar b ; \bar a \bar b
      if(i==j && a/L == b/L && a!=b)
         result += mol.getV(a_,b_,c_,d_) - mol.getV(a_,b_,d_,c_);

      // a \bar b ; a \bar b
      // \bar a b ; \bar a b
      if(i==j && a/L != b/L && a%L!=b%L)
         result += mol.getV(a_,b_, c_,d_);

      return result;
   };

   for(int i=0;i<block->getn();++i)
      for(int j=i;j<block->getn();++j)
         (*block)(i,j) = (*block)(j,i) = calc_elem(i,j);

   for(int i=0;i<diag.size();i++)
      // keep in mind that the degen of the vector is 4. We need prefactor of 2, so
      // we end up with 0.5
      diag[i] = 0.5*calc_elem(L+i,L+i) + 0.5*calc_elem(L*L+i,L*L+i);
}

/**
 * Calculate the dot product with this DM2 and other one.
 * @param x the other DM2
 * @return the dot product between this and x
 */
double DM2::Dot(const DM2 &x) const
{
   double result = 0;

   int n = block->getn()*block->getn();
   int inc = 1;

   result += ddot_(&n,block->getpointer(),&inc,x.block->getpointer(),&inc);

   n = diag.size();
   result += 4 * ddot_(&n,diag.data(),&inc,x.diag.data(),&inc);

   return result;
}

/**
 * @return the trace of the this DM2 object
 */
double DM2::Trace() const
{
   double result = 0;
   for(auto elem : diag)
      result += elem;

   return block->trace() + 4 * result;
}

/**
 * Find the minimum with Newton-Raphson for the angle of a jacobi rotation between
 * orbitals k and l in the DOCI space.
 * @param k the first orbital
 * @param l the second orbital
 * @param start_angle the starting point for the Newton-Raphson (defaults to zero)
 * @param T function that returns the one-particle matrix elements
 * @param V function that returns the two-particle matrix elements
 * @return pair of the angle with the lowest energy and boolean, true => minimum, false => maximum
 */
std::pair<double,bool> DM2::find_min_angle(int k, int l, double start_angle, std::function<double(int,int)> &T, std::function<double(int,int,int,int)> &V) const
{
   assert(k!=l);

   const int L = block->getn();

   double theta = start_angle;

   const DM2 &rdm = *this;

   double cos2 = 2.0/(N-1.0)*(T(k,k)*rdm(k,k+L,k,k+L)+T(l,l)*rdm(l,l+L,l,l+L));

   double sin2 = 2.0/(N-1.0)*(T(l,l)*rdm(k,k+L,k,k+L)+T(k,k)*rdm(l,l+L,l,l+L));

   // 2sincos actually
   double sincos = 2.0/(N-1.0)*T(k,l)*(rdm(l,l+L,l,l+L)-rdm(k,k+L,k,k+L));

   for(int a=0;a<L;a++)
   {
      if(a==k || a==l)
         continue;

      cos2 += 2*V(k,k,a,a)*rdm(k,k+L,a,a+L)+2*V(l,l,a,a)*rdm(l,l+L,a,a+L)+2*(2*V(k,a,k,a)-V(k,a,a,k)+2.0/(N-1.0)*T(k,k))*rdm(k,a,k,a)+2*(2*V(l,a,l,a)-V(l,a,a,l)+2.0/(N-1.0)*T(l,l))*rdm(l,a,l,a);

      sin2 += 2*V(l,l,a,a)*rdm(k,k+L,a,a+L)+2*V(k,k,a,a)*rdm(l,l+L,a,a+L)+2*(2*V(k,a,k,a)-V(k,a,a,k)+2.0/(N-1.0)*T(k,k))*rdm(l,a,l,a)+2*(2*V(l,a,l,a)-V(l,a,a,l)+2.0/(N-1.0)*T(l,l))*rdm(k,a,k,a);

      sincos += 2*V(k,l,a,a)*(rdm(l,l+L,a,a+L)-rdm(k,k+L,a,a+L))+2*(2*V(k,a,l,a)-V(k,a,a,l)+2.0/(N-1.0)*T(k,l))*(rdm(l,a,l,a)-rdm(k,a,k,a));
   }

   const double cos4 = V(k,k,k,k)*rdm(k,k+L,k,k+L)+V(l,l,l,l)*rdm(l,l+L,l,l+L)+2*V(k,k,l,l)*rdm(k,k+L,l,l+L)+2*(2*V(k,l,k,l)-V(k,k,l,l))*rdm(k,l,k,l);

   const double sin4 = V(k,k,k,k)*rdm(l,l+L,l,l+L)+V(l,l,l,l)*rdm(k,k+L,k,k+L)+2*V(k,k,l,l)*rdm(k,k+L,l,l+L)+2*(2*V(k,l,k,l)-V(k,k,l,l))*rdm(k,l,k,l);

   // 2 x
   const double cos2sin2 = (2*V(k,k,l,l)+V(k,l,k,l))*(rdm(k,k+L,k,k+L)+rdm(l,l+L,l,l+L))+((V(k,k,k,k)+V(l,l,l,l)-2*(V(k,l,k,l)+V(k,k,l,l))))*rdm(k,k+L,l,l+L)+(V(k,k,k,k)+V(l,l,l,l)-6*V(k,k,l,l)+2*V(k,l,k,l))*rdm(k,l,k,l);

   // 4 x
   const double sin3cos = V(k,l,k,k)*rdm(l,l+L,l,l+L)-V(k,l,l,l)*rdm(k,k+L,k,k+L)-(V(k,l,k,k)-V(k,l,l,l))*(rdm(k,k+L,l,l+L)+rdm(k,l,k,l));

   // 4 x
   const double cos3sin = V(k,l,l,l)*rdm(l,l+L,l,l+L)-V(k,l,k,k)*rdm(k,k+L,k,k+L)+(V(k,l,k,k)-V(k,l,l,l))*(rdm(k,k+L,l,l+L)+rdm(k,l,k,l));

   // A*cos(t)^4+B*sin(t)^4+C*cos(t)^2+D*sin(t)^2+2*E*cos(t)*sin(t)+2*F*cos(t)^2*sin(t)^2+4*G*sin(t)*cos(t)^3+4*H*sin(t)^3*cos(t)

   // (16*G-16*H)*cos(t)^4+(-4*A-4*B+8*F)*sin(t)*cos(t)^3+(4*E-12*G+20*H)*cos(t)^2+(4*B-2*C+2*D-4*F)*sin(t)*cos(t)-2*E-4*H
   auto gradient = [&] (double theta) -> double { 
      double cos = std::cos(theta);
      double sin = std::sin(theta);

      double res = 16*(cos3sin-sin3cos)*cos*cos*cos*cos - 4*(cos4+sin4-2*cos2sin2)*sin*cos*cos*cos + 4*(sincos-3*cos3sin+5*sin3cos)*cos*cos + 2*(2*sin4-cos2+sin2-2*cos2sin2)*sin*cos-2*sincos-4*sin3cos;

      return res;
   };

   // (-16*A-16*B+32*F)*cos(t)^4+(-64*G+64*H)*sin(t)*cos(t)^3+(12*A+20*B-4*C+4*D-32*F)*cos(t)^2+(-8*E+24*G-40*H)*sin(t)*cos(t)-4*B+2*C-2*D+4*F
   auto hessian = [&] (double theta) -> double { 
      double cos = std::cos(theta);
      double sin = std::sin(theta);

      double res = -16*(cos4+sin4-2*cos2sin2)*cos*cos*cos*cos+64*(sin3cos-cos3sin)*sin*cos*cos*cos+4*(3*cos4+5*sin4-cos2+sin2-8*cos2sin2)*cos*cos-8*(sincos-3*cos3sin+5*sin3cos)*sin*cos+2*(cos2-2*sin4-sin2+2*cos2sin2);
      
      return res;
   };

//   int Na = 5000;
//   for(int i=0;i<Na;i++)
//   {
//      double t = 2.0 * M_PI / (1.0*Na) * i;
//      std::cout << t << "\t" << calc_rotate(k,l,t,T,V)  << "\t" << gradient(t) << "\t" << hessian(t) << std::endl;
//   }

   const int max_iters = 20;
   const double convergence = 1e-12;

   double change = gradient(theta)*theta+hessian(theta)*theta*theta/2.0;

   // if it goes uphill, try flipping sign and try again
   if(change>0)
      theta *= -1;

   int iter;
   for(iter=0;iter<max_iters;iter++)
   {
      double dx = gradient(theta)/hessian(theta);

      theta -= dx;

      if(fabs(dx) < convergence)
         break;
   }

   if(iter>=max_iters)
      std::cout << "Reached max iters in find_min_angle: " << k << " " << l << std::endl;

   if(hessian(theta)<0)
      std::cout << "Found max!" << std::endl;

   return std::make_pair(theta, hessian(theta)>0);
}

/**
 * Calculate the energy change when you rotate orbital k and l over an angle of theta with
 * the rotation in the full space of the orbitals.
 * @param k the first orbital
 * @param l the second orbital
 * @param theta the angle to rotate over
 * @param T function that returns the one-particle matrix elements
 * @param V function that returns the two-particle matrix elements
 * @return the new energy
 */
double DM2::calc_rotate(int k, int l, double theta, std::function<double(int,int)> &T, std::function<double(int,int,int,int)> &V) const
{
   assert(k!=l);

   const int L = block->getn();

   const DM2 &rdm = *this;

   double energy = 4/(N-1.0)*(T(k,k)+T(l,l)) * rdm(k,l,k,l);

   double cos2 = 2.0/(N-1.0)*(T(k,k)*rdm(k,k+L,k,k+L)+T(l,l)*rdm(l,l+L,l,l+L));

   double sin2 = 2.0/(N-1.0)*(T(l,l)*rdm(k,k+L,k,k+L)+T(k,k)*rdm(l,l+L,l,l+L));

   // 2sincos actually
   double sincos = 2.0/(N-1.0)*T(k,l)*(rdm(l,l+L,l,l+L)-rdm(k,k+L,k,k+L));

   for(int a=0;a<L;a++)
   {
      if(a==k || a==l)
         continue;

      energy += 2.0/(N-1.0) * T(a,a) * (rdm(a,a+L,a,a+L)+2*rdm(a,k,a,k)+2*rdm(a,l,a,l));

      for(int b=0;b<L;b++)
      {
         if(b==k || b==l)
            continue;

         energy += 2.0/(N-1.0) * (T(a,a)+T(b,b)) * rdm(a,b,a,b);

         energy += V(a,a,b,b) * rdm(a,a+L,b,b+L);

         energy += (2*V(a,b,a,b)-V(a,b,b,a)) * rdm(a,b,a,b);
      }

      cos2 += 2*V(k,k,a,a)*rdm(k,k+L,a,a+L)+2*V(l,l,a,a)*rdm(l,l+L,a,a+L)+2*(2*V(k,a,k,a)-V(k,a,a,k)+2.0/(N-1.0)*T(k,k))*rdm(k,a,k,a)+2*(2*V(l,a,l,a)-V(l,a,a,l)+2.0/(N-1.0)*T(l,l))*rdm(l,a,l,a);

      sin2 += 2*V(l,l,a,a)*rdm(k,k+L,a,a+L)+2*V(k,k,a,a)*rdm(l,l+L,a,a+L)+2*(2*V(k,a,k,a)-V(k,a,a,k)+2.0/(N-1.0)*T(k,k))*rdm(l,a,l,a)+2*(2*V(l,a,l,a)-V(l,a,a,l)+2.0/(N-1.0)*T(l,l))*rdm(k,a,k,a);

      sincos += 2*V(k,l,a,a)*(rdm(l,l+L,a,a+L)-rdm(k,k+L,a,a+L))+2*(2*V(k,a,l,a)-V(k,a,a,l)+2.0/(N-1.0)*T(k,l))*(rdm(l,a,l,a)-rdm(k,a,k,a));
   }

   const double cos4 = V(k,k,k,k)*rdm(k,k+L,k,k+L)+V(l,l,l,l)*rdm(l,l+L,l,l+L)+2*V(k,k,l,l)*rdm(k,k+L,l,l+L)+2*(2*V(k,l,k,l)-V(k,k,l,l))*rdm(k,l,k,l);

   const double sin4 = V(k,k,k,k)*rdm(l,l+L,l,l+L)+V(l,l,l,l)*rdm(k,k+L,k,k+L)+2*V(k,k,l,l)*rdm(k,k+L,l,l+L)+2*(2*V(k,l,k,l)-V(k,k,l,l))*rdm(k,l,k,l);

   // 2 x
   const double cos2sin2 = (2*V(k,k,l,l)+V(k,l,k,l))*(rdm(k,k+L,k,k+L)+rdm(l,l+L,l,l+L))+((V(k,k,k,k)+V(l,l,l,l)-2*(V(k,l,k,l)+V(k,k,l,l))))*rdm(k,k+L,l,l+L)+(V(k,k,k,k)+V(l,l,l,l)-6*V(k,k,l,l)+2*V(k,l,k,l))*rdm(k,l,k,l);

   // 4 x
   const double sin3cos = V(k,l,k,k)*rdm(l,l+L,l,l+L)-V(k,l,l,l)*rdm(k,k+L,k,k+L)-(V(k,l,k,k)-V(k,l,l,l))*(rdm(k,k+L,l,l+L)+rdm(k,l,k,l));

   // 4 x
   const double cos3sin = V(k,l,l,l)*rdm(l,l+L,l,l+L)-V(k,l,k,k)*rdm(k,k+L,k,k+L)+(V(k,l,k,k)-V(k,l,l,l))*(rdm(k,k+L,l,l+L)+rdm(k,l,k,l));

   const double cos = std::cos(theta);
   const double sin = std::sin(theta);

   energy += cos*cos*cos*cos*cos4;

   energy += sin*sin*sin*sin*sin4;

   energy += cos*cos*cos2;

   energy += sin*sin*sin2;

   energy += 2*sin*cos*sincos;

   energy += 2*cos*cos*sin*sin*cos2sin2;

   energy += 4*cos*sin*sin*sin*sin3cos;

   energy += 4*cos*cos*cos*sin*cos3sin;

   return energy;
}

/**
 * Calculate S^2
 * @return S(S+1) of the current RDM
 */
double DM2::S2() const
{
   double S2 = 0;
   auto n_sp = block->getn();
   auto L = n_sp;
   auto M = 2*n_sp;
   auto n_tp = M*(M-1)/2;

    for(int i=0;i<n_tp;++i)
    {
        int a = (*tp2sp)(i,0);
        int b = (*tp2sp)(i,1);

        double s_a = ( 1.0 - 2 * (a / L) )/2;
        double s_b = ( 1.0 - 2 * (b / L) )/2;

        S2 += ( (1.0 + s_a*s_a + s_b*s_b)/(N - 1.0) + 2*s_a*s_b ) * (*this)(a,b,a,b);
    }

   for(int a=0;a<L;a++)
      for(int b=0;b<L;b++)
         S2 += (*this)(a,L+b,a+L,b);

   return S2;
}

/**
 * Calculate Sz
 * @return Sz of the current RDM
 */
double DM2::Sz() const
{
   double Sz = 0;
   auto n_sp = block->getn();
   auto L = n_sp;
   auto M = 2*n_sp;
   auto n_tp = M*(M-1)/2;

    for(int i=0;i<n_tp;++i)
    {
        int a = (*tp2sp)(i,0);
        int b = (*tp2sp)(i,1);

        double s_a = ( 1.0 - 2 * (a / L) )/2;
        double s_b = ( 1.0 - 2 * (b / L) )/2;

        Sz += 1.0/(N-1.0) * (s_a + s_b) * (*this)(a,b,a,b);
    }

   return Sz;
}

/* vim: set ts=3 sw=3 expandtab :*/
