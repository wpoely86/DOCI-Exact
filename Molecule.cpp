#include <iostream>
#include <assert.h>
#include <hdf5.h>

#include "Molecule.h"

// this helps to check the return codes of HDF5 calls
#define HDF5_STATUS_CHECK(status) if(status < 0) std::cerr << __FILE__ << ":" << __LINE__ << ": Problem with writing to file. Status code=" << status << std::endl;

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
Molecule::Molecule(std::string filename)
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
Molecule::Molecule(const Molecule &orig)
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
Molecule::Molecule(Molecule &&orig)
{
   OEI = std::move(orig.OEI);
   TEI = std::move(orig.TEI);

   n_electrons = orig.n_electrons;
   nucl_rep = orig.nucl_rep;
   n_sp = orig.n_sp;
}

Molecule& Molecule::operator=(const Molecule &orig)
{
   OEI.reset(new helpers::matrix(*orig.OEI));
   TEI.reset(new helpers::matrix(*orig.TEI));

   n_electrons = orig.n_electrons;
   nucl_rep = orig.nucl_rep;
   n_sp = orig.n_sp;

   return *this;
}

Molecule& Molecule::operator=(Molecule &&orig)
{
   OEI = std::move(orig.OEI);
   TEI = std::move(orig.TEI);

   n_electrons = orig.n_electrons;
   nucl_rep = orig.nucl_rep;
   n_sp = orig.n_sp;

   return *this;
}

/**
 * Get the matrix element <a|T|b> where T is the
 * one body operator
 * @param a the first sp index
 * @param b the second sp index
 * @return <a|T|b>
 */
double Molecule::getT(int a, int b) const
{
   assert(a<n_sp && b<n_sp);

   return (*OEI)(a,b);
}

/**
 * Get the matrix element <ab|V|cd> where V is the
 * two body operator
 * @param a the first sp index
 * @param b the second sp index
 * @param c the thirth sp index
 * @param d the fourth sp index
 * @return <ab|V|cd>
 */
double Molecule::getV(int a, int b, int c, int d) const
{
   assert(a<n_sp && b<n_sp && c<n_sp && d<n_sp);

   return (*TEI)(a*n_sp + b, c*n_sp + d);
}


/**
 * @return the nuclear replusion energy
 */
double Molecule::get_nucl_rep() const
{
   return nucl_rep;
}

/**
 * @return size of the single particle basis
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

/* vim: set ts=3 sw=3 expandtab :*/
