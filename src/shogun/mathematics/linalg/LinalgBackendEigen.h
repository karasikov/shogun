/*
 * Copyright (c) 2016, Shogun-Toolbox e.V. <shogun-team@shogun-toolbox.org>
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the copyright holder nor the names of its
 *     contributors may be used to endorse or promote products derived from
 *     this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * Authors: 2016 Pan Deng, Soumyajit De, Heiko Strathmann, Viktor Gal
 */

#include <shogun/lib/SGVector.h>
#include <shogun/mathematics/eigen3.h>
#include <shogun/mathematics/linalg/LinalgBackendBase.h>

#ifndef LINALG_BACKEND_EIGEN_H__
#define LINALG_BACKEND_EIGEN_H__

namespace shogun
{

/** @brief linalg methods with Eigen3 backend */
class LinalgBackendEigen : public LinalgBackendBase
{
public:
	#define DEFINE_FOR_ALL_PTYPE(METHODNAME, Container) \
	METHODNAME(bool, Container); \
	METHODNAME(char, Container); \
	METHODNAME(int8_t, Container); \
	METHODNAME(uint8_t, Container); \
	METHODNAME(int16_t, Container); \
	METHODNAME(uint16_t, Container); \
	METHODNAME(int32_t, Container); \
	METHODNAME(uint32_t, Container); \
	METHODNAME(int64_t, Container); \
	METHODNAME(uint64_t, Container); \
	METHODNAME(float32_t, Container); \
	METHODNAME(float64_t, Container); \
	METHODNAME(floatmax_t, Container); \
	METHODNAME(complex128_t, Container); \

	#define DEFINE_FOR_REAL_PTYPE(METHODNAME, Container) \
	METHODNAME(bool, Container); \
	METHODNAME(char, Container); \
	METHODNAME(int8_t, Container); \
	METHODNAME(uint8_t, Container); \
	METHODNAME(int16_t, Container); \
	METHODNAME(uint16_t, Container); \
	METHODNAME(int32_t, Container); \
	METHODNAME(uint32_t, Container); \
	METHODNAME(int64_t, Container); \
	METHODNAME(uint64_t, Container); \
	METHODNAME(float32_t, Container); \
	METHODNAME(float64_t, Container); \
	METHODNAME(floatmax_t, Container);

	/** Implementation of @see LinalgBackendBase::add */
	#define BACKEND_GENERIC_ADD(Type, Container) \
	virtual Container<Type> add(const Container<Type>& a, const Container<Type>& b, Type alpha, Type beta) const \
	{  \
		return add_impl(a, b, alpha, beta); \
	}
	DEFINE_FOR_ALL_PTYPE(BACKEND_GENERIC_ADD, SGVector)
	#undef BACKEND_GENERIC_ADD

	/** Implementation of @see LinalgBackendBase::dot */
	#define BACKEND_GENERIC_DOT(Type, Container) \
	virtual Type dot(const Container<Type>& a, const Container<Type>& b) const \
	{  \
		return dot_impl(a, b);  \
	}
	DEFINE_FOR_ALL_PTYPE(BACKEND_GENERIC_DOT, SGVector)
	#undef BACKEND_GENERIC_DOT

	/** Implementation of @see LinalgBackendBase::mean */
	#define BACKEND_GENERIC_REAL_MEAN(Type, Container) \
	virtual float64_t mean(const Container<Type>& a) const \
	{  \
		return mean_impl(a);  \
	}
	DEFINE_FOR_REAL_PTYPE(BACKEND_GENERIC_REAL_MEAN, SGVector)
	DEFINE_FOR_REAL_PTYPE(BACKEND_GENERIC_REAL_MEAN, SGMatrix)
	#undef BACKEND_GENERIC_REAL_MEAN

	/** Implementation of @see LinalgBackendBase::mean */
	#define BACKEND_GENERIC_COMPLEX_MEAN(Container) \
	virtual complex128_t mean(const Container<complex128_t>& a) const \
	{  \
		return mean_impl(a);  \
	}
	BACKEND_GENERIC_COMPLEX_MEAN(SGVector)
	BACKEND_GENERIC_COMPLEX_MEAN(SGMatrix)
	#undef BACKEND_GENERIC_COMPLEX_MEAN

	/** Implementation of @see LinalgBackendBase::sum */
	#define BACKEND_GENERIC_SUM(Type, Container) \
	virtual Type sum(const Container<Type>& a, bool no_diag) const \
	{  \
		return sum_impl(a, no_diag);  \
	}
	DEFINE_FOR_ALL_PTYPE(BACKEND_GENERIC_SUM, SGVector)
	DEFINE_FOR_ALL_PTYPE(BACKEND_GENERIC_SUM, SGMatrix)
	#undef BACKEND_GENERIC_SUM

	/** Implementation of @see LinalgBackendBase::colwise_sum */
	#define BACKEND_GENERIC_COLWISE_SUM(Type, Container) \
	virtual SGVector<Type> colwise_sum(const Container<Type>& a, bool no_diag) const \
	{  \
		return colwise_sum_impl(a, no_diag); \
	}
	DEFINE_FOR_ALL_PTYPE(BACKEND_GENERIC_COLWISE_SUM, SGMatrix)
	#undef BACKEND_GENERIC_COLWISE_SUM

	/** Implementation of @see LinalgBackendBase::rowwise_sum */
	#define BACKEND_GENERIC_ROWWISE_SUM(Type, Container) \
	virtual SGVector<Type> rowwise_sum(const Container<Type>& a, bool no_diag) const \
	{  \
		return rowwise_sum_impl(a, no_diag); \
	}
	DEFINE_FOR_ALL_PTYPE(BACKEND_GENERIC_ROWWISE_SUM, SGMatrix)
	#undef BACKEND_GENERIC_ROWWISE_SUM

	#undef DEFINE_FOR_ALL_PTYPE

private:
	/** Eigen3 vector C = alpha*A + beta*B method */
	template <typename T>
	SGVector<T> add_impl(const SGVector<T>& a, const SGVector<T>& b, T alpha, T beta) const
	{
		SGVector<T> c(a.vlen);
		typename SGVector<T>::EigenVectorXtMap a_eig = a;
		typename SGVector<T>::EigenVectorXtMap b_eig = b;
		typename SGVector<T>::EigenVectorXtMap c_eig = c;

		c_eig = alpha * a_eig + beta * b_eig;
		return c;
	}

	/** Eigen3 vector dot-product method */
	template <typename T>
	T dot_impl(const SGVector<T>& a, const SGVector<T>& b) const
	{
		return (typename SGVector<T>::EigenVectorXtMap(a)).dot(typename SGVector<T>::EigenVectorXtMap(b));
	}

	/** Real eigen3 vector and matrix mean method */
	template <typename T, template <typename> class Container>
	typename std::enable_if<!std::is_same<T, complex128_t>::value, float64_t>::type
	mean_impl(const Container<T>& a) const
	{
		return sum_impl(a)/(float64_t(a.size()));
	}

	/** Complex eigen3 vector and matrix mean method */
	template<template <typename> class Container>
	complex128_t mean_impl(const Container<complex128_t>& a) const
	{
		return sum_impl(a)/(complex128_t(a.size()));
	}

	/** Eigen3 vector sum method */
	template <typename T>
	T sum_impl(const SGVector<T>& vec, bool no_diag=false) const
	{
		return (typename SGVector<T>::EigenVectorXtMap(vec)).sum();
	}

	/** Eigen3 matrix sum method */
	template <typename T>
	T sum_impl(const SGMatrix<T>& mat, bool no_diag=false) const
	{
		typename SGMatrix<T>::EigenMatrixXtMap m = mat;
		T sum = m.sum();
		if (no_diag)
			sum -= m.diagonal().sum();

		return sum;
	}

	/** Eigen3 matrix colwise sum method */
	template <typename T>
	SGVector<T> colwise_sum_impl(const SGMatrix<T>& mat, bool no_diag) const
	{
		SGVector<T> result(mat.num_cols);

		typename SGMatrix<T>::EigenMatrixXtMap mat_eig = mat;
		typename SGVector<T>::EigenVectorXtMap result_eig = result;

		result_eig = mat_eig.colwise().sum();

		// remove the main diagonal elements if required
		if (no_diag)
		{
			index_t len_major_diag = mat_eig.rows() < mat_eig.cols()
				? mat_eig.rows() : mat_eig.cols();
			for (index_t i = 0; i < len_major_diag; ++i)
				result_eig[i] -= mat_eig(i,i);
		}

		return result;
	}

	/** Eigen3 matrix rowwise sum method */
	template <typename T>
	SGVector<T> rowwise_sum_impl(const SGMatrix<T>& mat, bool no_diag) const
	{
		SGVector<T> result(mat.num_rows);

		typename SGMatrix<T>::EigenMatrixXtMap mat_eig = mat;
		typename SGVector<T>::EigenVectorXtMap result_eig = result;

		result_eig = mat_eig.rowwise().sum();

		// remove the main diagonal elements if required
		if (no_diag)
		{
			index_t len_major_diag = mat_eig.rows() < mat_eig.cols()
				? mat_eig.rows() : mat_eig.cols();
			for (index_t i = 0; i < len_major_diag; ++i)
				result_eig[i] -= mat_eig(i,i);
		}

		return result;
	}

};

}

#endif //LINALG_BACKEND_EIGEN_H__
