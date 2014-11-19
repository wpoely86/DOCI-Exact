#include <iostream>
#include <assert.h>
#include <hdf5.h>

#include "Molecule.h"
#include "Permutation.h"

using namespace doci;

// this helps to check the return codes of HDF5 calls
#define HDF5_STATUS_CHECK(status) if(status < 0) std::cerr << __FILE__ << ":" << __LINE__ << ": Problem with writing to file. Status code=" << status << std::endl;

/**
 * @return the nuclear replusion energy
 */
double Molecule::get_nucl_rep() const
{
   return nucl_rep;
}

/**
 * @return size of the single particle basis (without spin)
 */
unsigned int Molecule::get_n_sp() const
{
   return n_sp;
}

/**
 * @return number of electrons
 */
unsigned int Molecule::get_n_electrons() const
{
   return n_electrons;
}



/**
 * Constructor. This reads a HDF5 file with the integrals. It expects following
 * file format:
 * /integrals/OEI => array with One electron integrals
 * /integrals/TEI => array with Two electron integrals
 * /integrals should have following attributes:
 * - nelectrons => number of electrons
 * - nuclear_repulsion_energy => nuclear repulsion energy
 * - sp_dim => size of the single particles basis
 * @param filename the HDF5 file to read the OEI and TEI from
 */
PSI_C1_Molecule::PSI_C1_Molecule(std::string filename)
{
   hid_t       file_id, group_id, dataset_id, attribute_id;
   herr_t      status;

   file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   HDF5_STATUS_CHECK(file_id);

   group_id = H5Gopen(file_id, "/integrals", H5P_DEFAULT);
   HDF5_STATUS_CHECK(group_id);

   attribute_id = H5Aopen(group_id, "nelectrons", H5P_DEFAULT);
   HDF5_STATUS_CHECK(attribute_id);

   status = H5Aread(attribute_id, H5T_NATIVE_UINT, &n_electrons);
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   attribute_id = H5Aopen(group_id, "sp_dim", H5P_DEFAULT);
   HDF5_STATUS_CHECK(attribute_id);

   status = H5Aread(attribute_id, H5T_NATIVE_UINT, &n_sp);
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   if(n_sp > Permutation::getMax() )
      throw std::overflow_error("The used type is not big enough to store all single particle states!");

   attribute_id = H5Aopen(group_id, "nuclear_repulsion_energy", H5P_DEFAULT);
   HDF5_STATUS_CHECK(attribute_id);

   status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &nucl_rep);
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   // allocate space
   OEI.reset(new helpers::matrix(n_sp, n_sp));

   TEI.reset(new helpers::matrix(n_sp*n_sp, n_sp*n_sp));

   dataset_id = H5Dopen(group_id, "OEI", H5P_DEFAULT);
   HDF5_STATUS_CHECK(dataset_id);

   status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, OEI->getpointer());
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);

   dataset_id = H5Dopen(group_id, "TEI", H5P_DEFAULT);
   HDF5_STATUS_CHECK(dataset_id);

   status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, TEI->getpointer());
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);

   status = H5Gclose(group_id);
   HDF5_STATUS_CHECK(status);

   status = H5Fclose(file_id);
   HDF5_STATUS_CHECK(status);
}

/**
 * Copy constructor
 */
PSI_C1_Molecule::PSI_C1_Molecule(const PSI_C1_Molecule &orig)
{
   OEI.reset(new helpers::matrix(*orig.OEI));
   TEI.reset(new helpers::matrix(*orig.TEI));

   n_electrons = orig.n_electrons;
   nucl_rep = orig.nucl_rep;
   n_sp = orig.n_sp;
}

/**
 * Move constructor
 */
PSI_C1_Molecule::PSI_C1_Molecule(PSI_C1_Molecule &&orig)
{
   OEI = std::move(orig.OEI);
   TEI = std::move(orig.TEI);

   n_electrons = orig.n_electrons;
   nucl_rep = orig.nucl_rep;
   n_sp = orig.n_sp;
}

PSI_C1_Molecule& PSI_C1_Molecule::operator=(const PSI_C1_Molecule &orig)
{
   OEI.reset(new helpers::matrix(*orig.OEI));
   TEI.reset(new helpers::matrix(*orig.TEI));

   n_electrons = orig.n_electrons;
   nucl_rep = orig.nucl_rep;
   n_sp = orig.n_sp;

   return *this;
}

PSI_C1_Molecule& PSI_C1_Molecule::operator=(PSI_C1_Molecule &&orig)
{
   OEI = std::move(orig.OEI);
   TEI = std::move(orig.TEI);

   n_electrons = orig.n_electrons;
   nucl_rep = orig.nucl_rep;
   n_sp = orig.n_sp;

   return *this;
}

/**
 * Clone the current object. Needed because
 * the base class is pure virtual and thus has
 * no copy constructor.
 * @return a clone of this object
 */
PSI_C1_Molecule* PSI_C1_Molecule::clone() const
{
   return new PSI_C1_Molecule(*this);
}

/**
 * Get the matrix element \f$<a|\hat T|b>\f$ where T is the
 * one body operator
 * @param a the first sp index
 * @param b the second sp index
 * @return \f$<a|\hat T|b>\f$
 */
double PSI_C1_Molecule::getT(int a, int b) const
{
   assert(a<n_sp && b<n_sp);

   return (*OEI)(a,b);
}

/**
 * Get the matrix element \f$<ab|\hat V|cd>\f$ where V is the
 * two body operator
 * @param a the first sp index
 * @param b the second sp index
 * @param c the thirth sp index
 * @param d the fourth sp index
 * @return \f$<ab|\hat V|cd>\f$
 */
double PSI_C1_Molecule::getV(int a, int b, int c, int d) const
{
   assert(a<n_sp && b<n_sp && c<n_sp && d<n_sp);

   return (*TEI)(a*n_sp + b, c*n_sp + d);
}

/**
 * Print the molecular integrals. Usefull for input into
 * other codes
 */
void PSI_C1_Molecule::Print() const
{
   auto L = get_n_sp();

   printf("%20.15f\t%d\t%d\t0\t0\n", 0.0, 0, 0);
   for(int a=0;a<L;a++)
      for(int b=0;b<L;b++)
         printf("%20.15f\t%d\t%d\t0\t0\n", (*OEI)(a,b), a+1,b+1);

   for(int a=0;a<L;a++)
      for(int b=0;b<L;b++)
         for(int c=0;c<L;c++)
            for(int d=0;d<L;d++)
               printf("%20.15f\t%d\t%d\t%d\t%d\n", (*TEI)(a*L+b,c*L+d), a+1,c+1,b+1,d+1);
}

/**
 * This will calculate the (restricted) Hartree-Fock energy.
 * @warning This will only work if the molecular orbitals are sorted
 * according to energy! Not if they are sorted according to irrep (as 
 * PSI does by default).
 * @return the (restricted) Hartree-Fock energy
 */
double PSI_C1_Molecule::HF_Energy() const
{
   double energy = 0;

   // one particle terms
   for(int a=0;a<n_electrons/2;a++)
      energy += 2 * getT(a,a);

   // two particle terms
   for(int a=0;a<n_electrons/2;a++)
      for(int b=0;b<n_electrons/2;b++)
         energy += 2 * getV(a,b,a,b) - getV(a,b,b,a);

   return energy;
}

/* vim: set ts=3 sw=3 expandtab :*/
