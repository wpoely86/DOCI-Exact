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

#include <iostream>
#include <algorithm>
#include <sstream>
#include <assert.h>
#include <hdf5.h>

#include "helpers.h"
#include "lapack.h"

// macro to help check return status of HDF5 functions
#define HDF5_STATUS_CHECK(status) {                 \
    if(status < 0)                                  \
    std::cerr << __FILE__ << ":" << __LINE__ <<     \
    ": Problem with writing to file. Status code="  \
    << status << std::endl;                         \
}

using namespace helpers;

/**
 * Empty matrix constructor. Don't use it
 * unless you know what you're doing
 */
matrix::matrix()
{
    this->n = 0;
    this->m = 0;
}

/**
 * @param n_ number of rows
 * @param m_ number of columns
 */
matrix::matrix(int n_, int m_)
{
    assert(n_ && m_);
    this->n = n_;
    this->m = m_;
    mat.reset(new double [n*m]);
}

/**
 * @param orig matrix to copy
 */
matrix::matrix(const matrix &orig)
{
    n = orig.n;
    m = orig.m;
    mat.reset(new double [n*m]);
    std::memcpy(mat.get(), orig.getpointer(), n*m*sizeof(double));
}

/**
 * move constructor
 * @param orig matrix to copy (descrutive)
 */
matrix::matrix(matrix &&orig)
{
    n = orig.n;
    m = orig.m;
    mat = std::move(orig.mat);
}

matrix& matrix::operator=(const matrix &orig)
{
    n = orig.n;
    m = orig.m;
    mat.reset(new double [n*m]);
    std::memcpy(mat.get(), orig.getpointer(), n*m*sizeof(double));
    return *this;
}

/**
 * Set all matrix elements equal to a value
 * @param val the value to use
 */
matrix& matrix::operator=(double val)
{
    for(int i=0;i<n*m;i++)
        mat[i] = val;

    return *this;
}

/**
 * @return number of rows
 */
int matrix::getn() const
{
    return n;
}

/**
 * @return number of columns
 */
int matrix::getm() const
{
    return m;
}

double matrix::operator()(int x,int y) const
{
    assert(x<n && y<m);
    return mat[x+y*n];
}

double& matrix::operator()(int x,int y)
{
    assert(x<n && y<m);
    return mat[x+y*n];
}

double& matrix::operator[](int x)
{
    assert(x<n*m);
    return mat[x];
}

double matrix::operator[](int x) const
{
    assert(x<n*m);
    return mat[x];
}

double* matrix::getpointer() const
{
    return mat.get();
}

/**
 * Matrix-Matrix product of A and B. Store result in this
 * @param A first matrix
 * @param B second matrix
 */
matrix& matrix::prod(matrix const &A, matrix const &B)
{
    char trans = 'N';

    double alpha = 1.0;
    double beta = 0.0;

    assert(A.n == n && B.m == m);

    dgemm_(&trans,&trans,&A.n,&B.m,&A.m,&alpha,A.mat.get(),&A.n,B.mat.get(),&B.n,&beta,mat.get(),&A.n);

    return *this;
}

/**
 * Do a SVD on this matrix and store left singular values in
 * this. Changes the size of the matrix!
 * @return list of singular values
 */
std::unique_ptr<double []> matrix::svd()
{
    char jobu = 'A';
    char jobvt = 'N';

    int count_sing = std::min(n,m);

    std::unique_ptr<double []> sing_vals(new double[count_sing]);

    // MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)).
    int lwork = std::max( 3*count_sing + std::max(n,m), 5*count_sing);
    std::unique_ptr<double []> work(new double[lwork]);

    std::unique_ptr<double []> vt(new double[n*n]);

    int info;

    dgesvd_(&jobu,&jobvt,&n,&m,mat.get(),&n,sing_vals.get(),vt.get(),&n,0,&m,work.get(),&lwork,&info);

    if(info)
        std::cerr << "svd failed. info = " << info << std::endl;

    // overwrite the matrix with the right singular vectors
    m = n;
    mat = std::move(vt);

    return sing_vals;
}

void matrix::Print() const
{
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            std::cout << i << " " << j << "\t" << (*this)(i,j) << std::endl;
}

double matrix::trace() const
{
    double result = 0;

    for(int i=0;i<std::min(m,n);i++)
        result += (*this)(i,i);

    return result;
}


void matrix::SaveToFile(std::string filename) const
{
    hid_t       file_id, dataset_id, dataspace_id, attribute_id;
    herr_t      status;

    file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_STATUS_CHECK(file_id);

    hsize_t dimarray = n*m;
    dataspace_id = H5Screate_simple(1, &dimarray, NULL);

    dataset_id = H5Dcreate(file_id, "matrix", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mat.get() );
    HDF5_STATUS_CHECK(status);

    status = H5Sclose(dataspace_id);
    HDF5_STATUS_CHECK(status);

    dataspace_id = H5Screate(H5S_SCALAR);

    attribute_id = H5Acreate (dataset_id, "n", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite (attribute_id, H5T_NATIVE_INT, &n );
    HDF5_STATUS_CHECK(status);

    status = H5Aclose(attribute_id);
    HDF5_STATUS_CHECK(status);
    
    attribute_id = H5Acreate (dataset_id, "m", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite (attribute_id, H5T_NATIVE_INT, &m );
    HDF5_STATUS_CHECK(status);

    status = H5Aclose(attribute_id);
    HDF5_STATUS_CHECK(status);

    status = H5Sclose(dataspace_id);
    HDF5_STATUS_CHECK(status);

    status = H5Dclose(dataset_id);
    HDF5_STATUS_CHECK(status);

    status = H5Fclose(file_id);
    HDF5_STATUS_CHECK(status);
}

void matrix::ReadFromFile(std::string filename) const
{
    hid_t       file_id, dataset_id, attribute_id;
    herr_t      status;
    int n_, m_;

    file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    HDF5_STATUS_CHECK(file_id);

    dataset_id = H5Dopen(file_id, "matrix", H5P_DEFAULT);
    HDF5_STATUS_CHECK(dataset_id);

    attribute_id = H5Aopen(dataset_id, "n", H5P_DEFAULT);
    HDF5_STATUS_CHECK(attribute_id);

    status = H5Aread(attribute_id, H5T_NATIVE_INT, &n_);
    HDF5_STATUS_CHECK(status);

    status = H5Aclose(attribute_id);
    HDF5_STATUS_CHECK(status);

    attribute_id = H5Aopen(dataset_id, "m", H5P_DEFAULT);
    HDF5_STATUS_CHECK(attribute_id);

    status = H5Aread(attribute_id, H5T_NATIVE_INT, &m_);
    HDF5_STATUS_CHECK(status);

    status = H5Aclose(attribute_id);
    HDF5_STATUS_CHECK(status);

    if(n!=n_ || m!=m_)
        std::cerr << "Matrix size not compatable with file: " << n_ << "x" << m_ << std::endl;
    else
    {
        status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mat.get());
        HDF5_STATUS_CHECK(status);
    }

    status = H5Dclose(dataset_id);
    HDF5_STATUS_CHECK(status);

    status = H5Fclose(file_id);
    HDF5_STATUS_CHECK(status);
}

/* vim: set ts=8 sw=4 tw=0 expandtab :*/
