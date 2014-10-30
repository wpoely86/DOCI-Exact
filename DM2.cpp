#include <hdf5.h>
#include <assert.h>

#include "DM2.h"
#include "Hamiltonian.h"

std::unique_ptr<helpers::matrix> DM2::sp2tp = nullptr;
std::unique_ptr<helpers::matrix> DM2::tp2sp = nullptr;

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

   this->n =  n;
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
   n = mol.get_n_electrons();
}

DM2::DM2(const DM2 &orig)
{
   assert(sp2tp && tp2sp);

   block.reset(new helpers::matrix(*orig.block));
   diag = orig.diag;
   n = orig.n;
}

DM2::DM2(DM2 &&orig)
{
   block = std::move(orig.block);
   diag = std::move(orig.diag);
   n = orig.n;
}

DM2& DM2::operator=(const DM2 &orig)
{
   block.reset(new helpers::matrix(*orig.block));
   diag = orig.diag;
   n = orig.n;

   return *this;
}

DM2& DM2::operator=(DM2 &&orig)
{
   block = std::move(orig.block);
   diag = std::move(orig.diag);
   n = orig.n;

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

   group_id = H5Gcreate(file_id, "/RDM", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   hsize_t dimblock = block->getn() * block->getn();

   dataspace_id = H5Screate_simple(1, &dimblock, NULL);

   dataset_id = H5Dcreate(group_id, "Block", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

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
   status = H5Awrite (attribute_id, H5T_NATIVE_UINT, &n );
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

   unsigned int n;
   status = H5Aread(attribute_id, H5T_NATIVE_UINT, &n);
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   DM2 dm2(L,n);

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
   return n;
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
   (*block) = 0;
   std::fill(diag.begin(), diag.end(), 0);

   perm.reset();
   auto& perm_bra = perm;

   for(unsigned int i=0;i<eigv.size();++i)
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
         (*block)(s,s) += eigv[i] * eigv[i];

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

            diag[idx] += eigv[i] * eigv[i];
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
            (*block)(r,s) += eigv[i] * eigv[j];
            (*block)(s,r) += eigv[i] * eigv[j];
         }
      }

      perm_bra.next();
   }

   auto tmp = block->trace();
   for(auto elem : diag)
      tmp += 4 * elem;
   std::cout << "Trace: " << tmp << std::endl;
   std::cout << "Trace 2: " << block->trace() << std::endl;
}

std::ostream &operator<<(std::ostream &output,DM2 &dm2)
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

/* vim: set ts=3 sw=3 expandtab :*/
