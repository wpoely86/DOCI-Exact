#include <cmath>
#include <hdf5.h>
#include "SparseMatrix_CRS.h"

// this helps to check the return codes of HDF5 calls
#define HDF5_STATUS_CHECK(status) if(status < 0) std::cerr << __FILE__ << ":" << __LINE__ << ": Problem with writing to file. Status code=" << status << std::endl;

using namespace helpers;

/**
 * Construct SparseMatrix_CRS object for n x n matrix
 * @param n the number of rows/columns
 */
SparseMatrix_CRS::SparseMatrix_CRS(unsigned int n)
{
    this->n = n;
    row.reserve(n+1);
}

/**
 * Read only access operator
 * @param i the row number
 * @param j the column number
 * @return the matrix element
 */
double SparseMatrix_CRS::operator()(unsigned int i,unsigned int j) const
{
   assert(i<n && j<n);

    for(unsigned int k=row[i];k<row[i+1];k++)
       if( col[k] == j )
          return data[k];

    return 0;
}

/**
 * @return the number of rows
 */
unsigned int SparseMatrix_CRS::gn() const
{
    return n;
}

/**
 * Convert a dense matrix to CRS format
 * @param dense the matrix to convert
 */
void SparseMatrix_CRS::ConvertFromMatrix(const helpers::matrix &dense)
{
   assert(dense.getm() == dense.getn());
   this->n = dense.getn();
   row.resize(n+1);

   data.clear();
   col.clear();

   row[0] = 0;

   for(unsigned int i=0;i<n;i++)
   {
      for(unsigned int j=0;j<n;j++)
         if( fabs(dense(i,j)) > 1e-14 )
         {
            data.push_back(dense(i,j));
            col.push_back(j);
         }

      row[i+1] = col.size();
   }

   row.back() = col.size();
}

/**
 * Convert this CRS matrix to a dense symmetric matrix (only works for square matrices)
 * @param dense the matrix to fill
 */
void SparseMatrix_CRS::ConvertToMatrix(helpers::matrix &dense) const
{
   assert(dense.getm() == dense.getn() && dense.getn() == n);
   dense = 0;

   for(unsigned int i=0;i<row.size()-1;i++)
      for(unsigned int k=row[i];k<row[i+1];k++)
         dense(i,col[k]) = dense(col[k],i) = data[k];
}

/**
 * Print the raw CRS data to stdout
 */
void SparseMatrix_CRS::PrintRaw() const
{
    std::cout << "Data(" << data.size() << "):" << std::endl;
    for(unsigned int i=0;i<data.size();i++)
        std::cout << data[i] << " ";
    std::cout << std::endl;

    std::cout << "Col indices:" << std::endl;
    for(unsigned int i=0;i<col.size();i++)
        std::cout << col[i] << " ";
    std::cout << std::endl;

    std::cout << "Row indices:" << std::endl;
    for(unsigned int i=0;i<row.size();i++)
        std::cout << row[i] << " ";
    std::cout << std::endl;
}

/**
 * Print sparse matrix to output
 * @param output the ostream to print to
 * @param matrix_p the matrix to print
 * @return the filled ostream (with the matrix)
 */
std::ostream &operator<<(std::ostream &output,helpers::SparseMatrix_CRS &matrix_p)
{
   for(unsigned int i=0;i<matrix_p.row.size()-1;i++)
      for(unsigned int k=matrix_p.row[i];k<matrix_p.row[i+1];k++)
         output << i << "\t" << matrix_p.col[k] << "\t" << matrix_p.data[k] << std::endl;

   return output;
}

/**
 * Adds a new column element to the current row.
 * To use this, first call NewRow() to start a row and then
 * use PushToRow() to add elements to that row. Always end
 * with calling NewRow() again.
 * @param j column
 * @param value the matrix element value
 */
void SparseMatrix_CRS::PushToRow(unsigned int j, double value)
{
   if(col.empty() || row.back() == col.size() || col.back() < j)
   {
      data.push_back(value);
      col.push_back(j);
   }
   else
   {
      unsigned int begin = row.back();
      for(unsigned int i=begin;i<col.size();i++)
      {
         if( col[i] > j )
         {
            col.insert(col.begin() + i,j);
            data.insert(data.begin() + i,value);
            break;
         } else if (col[i] == j)
         {
            data[i] += value;

            if(fabs(data[i])<1e-14)
            {
               data.erase(data.begin() + i);
               col.erase(col.begin() + i);
            }
            break;
         }
      }
   }
}

/**
 * Adds a new column element to the current row.
 * To use this, first call NewRow() to start a row and then
 * use PushToRowNext() to add elements to that row. Always end
 * with calling NewRow() again.
 *
 * IMPORTANT: the PushToRow is a general function, you can push
 * values to the row in any order. This method assumes you always
 * add a new element in a column next to the previous one.
 *
 * @param j column
 * @param value the matrix element value
 */
void SparseMatrix_CRS::PushToRowNext(unsigned int j, double value)
{
   assert(col.empty() || row.back() == col.size() || col.back() < j);

   data.push_back(value);
   col.push_back(j);
}

/**
 * Adds the next row to the sparsematrix
 */
void SparseMatrix_CRS::NewRow()
{
   if(row.size() == (n+1))
      return;

   row.push_back(data.size());

   // fill the lower part with data
   // from the upper part
   for(unsigned int j=1;j<row.size();j++)
      for(unsigned int k=row[j-1]+1;k<row[j];k++)
      {
         if(col[k] == (row.size()-1))
         {
            data.push_back(data[k]);
            col.push_back(j-1);
            break;
         }
      }
}

/**
 * Do the matrix vector product y = A * x + beta * y
 * @param x a m component vector
 * @param y a n component vector
 * @param beta the multiply factor for y
 */
void SparseMatrix_CRS::mvprod(const double *x, double *y, double beta) const
{
#pragma omp parallel
   for(unsigned int i=0;i<n;i++)
   {
      y[i] *= beta;

      for(unsigned int k=row[i];k<row[i+1];k++)
      {
         y[i] += data[k] * x[col[k]];
      }
   }
}

/**
 * Save a SparseMatrix_CRS to a HDF5 file
 * @param filename the name of the file to write to
 * @param name the name of the group in the HDF5 file
 * @param append add or overwrite file
 */
int SparseMatrix_CRS::WriteToFile(const char *filename, const char *name, bool append) const
{
   hid_t       file_id, group_id, dataset_id, attribute_id, dataspace_id, scalar_id;
   herr_t      status;

   if(append)
      file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
   else
      file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

   group_id = H5Gcreate(file_id, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   hsize_t dimblock = data.size();

   scalar_id = H5Screate(H5S_SCALAR);

   dataspace_id = H5Screate_simple(1, &dimblock, NULL);

   dataset_id = H5Dcreate(group_id, "data", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data() );
   HDF5_STATUS_CHECK(status);

   unsigned int size = data.size();
   attribute_id = H5Acreate (dataset_id, "size", H5T_STD_U64LE, scalar_id, H5P_DEFAULT, H5P_DEFAULT);
   HDF5_STATUS_CHECK(attribute_id);
   status = H5Awrite (attribute_id, H5T_NATIVE_UINT, &size );
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);

   dataset_id = H5Dcreate(group_id, "col", H5T_STD_U64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   status = H5Dwrite(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, col.data() );
   HDF5_STATUS_CHECK(status);

   size = col.size();
   attribute_id = H5Acreate (dataset_id, "size", H5T_STD_U64LE, scalar_id, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Awrite (attribute_id, H5T_NATIVE_UINT, &size );
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);

   status = H5Sclose(dataspace_id);
   HDF5_STATUS_CHECK(status);

   dimblock = row.size();

   dataspace_id = H5Screate_simple(1, &dimblock, NULL);

   dataset_id = H5Dcreate(group_id, "row", H5T_STD_U64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   status = H5Dwrite(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, row.data() );
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);

   status = H5Sclose(dataspace_id);
   HDF5_STATUS_CHECK(status);

   dataset_id = H5Dcreate(group_id, "n", H5T_STD_U64LE, scalar_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   status = H5Dwrite(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &n );
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);

   status = H5Sclose(scalar_id);
   HDF5_STATUS_CHECK(status);

   status = H5Gclose(group_id);
   HDF5_STATUS_CHECK(status);

   status = H5Fclose(file_id);
   HDF5_STATUS_CHECK(status);

   return 0;
}

/**
 * Read a SparseMatrix_CRS from a HDF5 file
 * @param filename the name of the file to read from
 * @param name the name of the group in the HDF5 file
 */
int SparseMatrix_CRS::ReadFromFile(const char *filename, const char *name)
{
   hid_t       file_id, group_id, dataset_id, attribute_id;
   herr_t      status;

   file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
   HDF5_STATUS_CHECK(file_id);

   group_id = H5Gopen(file_id, name, H5P_DEFAULT);
   HDF5_STATUS_CHECK(group_id);

   dataset_id = H5Dopen(group_id, "n", H5P_DEFAULT);
   HDF5_STATUS_CHECK(dataset_id);

   status = H5Dread(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &n);
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);

   row.resize(n+1);

   unsigned int size;

   dataset_id = H5Dopen(group_id, "data", H5P_DEFAULT);
   HDF5_STATUS_CHECK(dataset_id);

   attribute_id = H5Aopen(dataset_id, "size", H5P_DEFAULT);
   HDF5_STATUS_CHECK(attribute_id);

   status = H5Aread(attribute_id, H5T_NATIVE_UINT, &size);
   HDF5_STATUS_CHECK(status);


   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   data.resize(size);

   status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);


   dataset_id = H5Dopen(group_id, "col", H5P_DEFAULT);
   HDF5_STATUS_CHECK(dataset_id);

   attribute_id = H5Aopen(dataset_id, "size", H5P_DEFAULT);
   HDF5_STATUS_CHECK(attribute_id);

   status = H5Aread(attribute_id, H5T_NATIVE_UINT, &size);
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   col.resize(size);

   status = H5Dread(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, col.data());
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);


   dataset_id = H5Dopen(group_id, "row", H5P_DEFAULT);
   HDF5_STATUS_CHECK(dataset_id);

   status = H5Dread(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, row.data());
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);

   status = H5Gclose(group_id);
   HDF5_STATUS_CHECK(status);

   status = H5Fclose(file_id);
   HDF5_STATUS_CHECK(status);

   return 0;
}

/**
 * Get the number of elements in a row
 * @param idx the row to consider
 * @return the number of non-zero elements in row idx
 */
unsigned int SparseMatrix_CRS::NumOfElInRow(unsigned int idx) const
{
   return (row[idx+1]-row[idx]);
}

/**
 * Get a element out of a row
 * @param row_index the number of the row where to get an element of
 * @param element_index the index of the element (index of the non-zero elements, not the column index)
 * @return the value of the element
 */
double SparseMatrix_CRS::GetElementInRow(unsigned int row_index, unsigned int element_index) const
{
   return data[row[row_index]+element_index];
}

/**
 * Get the column number of a element in a row. Used together with GetElementInRow()
 * @param row_index the index of the row
 * @param element_index the index of the element in the row (index of the non-zero elements, not the column index)
 * @return the column number of the element
 */
unsigned int SparseMatrix_CRS::GetElementColIndexInRow(unsigned int row_index, unsigned int element_index) const
{
   return col[row[row_index]+element_index];
}


/**
 * Set a guess for the number of non-zero elements in
 * the matrix. This method will reserve enough memory to store
 * at least count elements
 * @param the number of guessed non-zero elements
 */
void SparseMatrix_CRS::SetGuess(unsigned int count)
{
   data.reserve(count);
   col.reserve(count);
}

/* vim: set ts=3 sw=3 expandtab :*/
