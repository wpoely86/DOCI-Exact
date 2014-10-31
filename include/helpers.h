/* Copyright (C) 2014  Ward Poelmans

   This file is part of Hubbard-GPU.

   Hubbard-GPU is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Hubbard-GPU is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Hubbard-GPU.  If not, see <http://www.gnu.org/licenses/>.
   */

#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <memory>
#include <vector>
#include <cstring>
#include <assert.h>
#include "SparseMatrix_CCS.h"

namespace helpers {

/**
 * Helper class, wrapper around a double array. Has methods to get the number of rows
 * and columns. A simple matrix class.
 */
class matrix
{
    public:
        matrix();

        matrix(int n_, int m_);

        matrix(const matrix &orig);

        matrix(matrix &&orig);

        virtual ~matrix() { }

        matrix& operator=(const matrix &orig);

        matrix& operator=(double val);

        int getn() const;

        int getm() const;

        double operator()(int x,int y) const;

        double& operator()(int x,int y);

        double& operator[](int x);

        double operator[](int x) const;

        double* getpointer() const;

        matrix& prod(matrix const &A, matrix const &B);

        std::unique_ptr<double []> svd();

        void mvprod(const double *x, double *y, double beta) const;

        void Print() const;

        double trace() const;

        std::vector<double> GetColumn(unsigned int) const;

        const double* GetColumnRaw(unsigned int idx) const;

        void SaveToFile(std::string filename) const;

        void ReadFromFile(std::string filename) const;

    private:
        //!n by m array of double
        std::unique_ptr<double []> mat;
        //! number of rows
        int n;
        //! number of columns
        int m;
};

}

#endif /* MATRIX_H */

/* vim: set ts=8 sw=4 tw=0 expandtab :*/
