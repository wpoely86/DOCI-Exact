#ifndef SPARSEMATRIX_CRS_H
#define SPARSEMATRIX_CRS_H

#include <iostream>
#include <vector>

#include "helpers.h"

// dark magic to get the friend operator<< to work...
namespace helpers { class SparseMatrix_CRS; };
std::ostream &operator<<(std::ostream &output,helpers::SparseMatrix_CRS &matrix_p);

namespace helpers {

/**
 * @author Ward Poelmans
 * @date 03-07-2012\n\n
 * This is a class written for sparse n x n matrices to use on the gpu. It uses the CRS format to store
 * a matrix. Only symmetric matrices!
 */

class SparseMatrix_CRS
{
   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << matrix << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << matrix << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param matrix_p de SparseMatrix_CRS you want to print
    */
   friend std::ostream &::operator<<(std::ostream &output,helpers::SparseMatrix_CRS &matrix_p);

   public:

      SparseMatrix_CRS(unsigned int n);

      virtual ~SparseMatrix_CRS() = default;

      //easy to access the numbers
      double operator()(unsigned int i,unsigned int j) const;

      unsigned int gn() const;

      void ConvertFromMatrix(const helpers::matrix &dense);

      void ConvertToMatrix(helpers::matrix &dense) const;

      void PrintRaw() const;

      void PushToRow(unsigned int j, double value);

      void PushToRowNext(unsigned int j, double value);

      void NewRow();

      void mvprod(const double *, double *) const;

      void mvprod(const double *, double *, double) const;

      void SetGuess(unsigned int);

      int WriteToFile(const char*,const char*,bool=false) const;

      int ReadFromFile(const char*,const char*);

      unsigned int NumOfElInRow(unsigned int idx) const;

      double GetElementInRow(unsigned int row_index, unsigned int element_index) const;

      unsigned int GetElementColIndexInRow(unsigned int row_index, unsigned int element_index) const;

   private:

      //! Array that holds the non zero values
      std::vector<double> data;
      //! Array that holds the column indexes
      std::vector<unsigned int> col;
      //! Array that holds the row index of data
      std::vector<unsigned int> row;

      //!dimension of the matrix (number of rows/columns)
      unsigned int n;
};

}

#endif /* SPARSEMATRIX_CRS_H */

/* vim: set ts=3 sw=3 expandtab :*/
