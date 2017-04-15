/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Written (W) 2013 Kevin Hughes
 *
 * Thanks to Andreas Ziehe and Cedric Gouy-Pailler
 */

#ifndef APPROXJOINTDIAGONALIZER_H_
#define APPROXJOINTDIAGONALIZER_H_

#include <shogun/lib/config.h>


#include <shogun/lib/common.h>
#include <shogun/lib/SGMatrix.h>
#include <shogun/lib/SGNDArray.h>
#include <shogun/base/SGObject.h>

#include <shogun/mathematics/Math.h>

namespace shogun
{

/** @brief Class ApproxJointDiagonalizer defines an
 * Approximate Joint Diagonalizer (AJD) interface.
 *
 * AJD finds the matrix V that best diagonalizes
 * a set \f${C^1 ... C^k}\f$ of real valued symmetric
 * \f$NxN\f$ matrices - \f$V*C*V^T\f$
 */
class CApproxJointDiagonalizer : public CSGObject
{
	public:

		/** constructor */
		CApproxJointDiagonalizer() : CSGObject()
		{
		};

		/** destructor */
		virtual ~CApproxJointDiagonalizer()
		{
		}

		/** Computes the matrix V that best diagonalizes C
		 * @param C the set of matrices to be diagonalized
		 * @param V0 an estimate of the matrix V
		 * @param eps machine epsilon or desired epsilon
		 * @param itermax maximum number of iterations
		 * @return V the matrix that best diagonalizes C
		 */
		virtual SGMatrix<float64_t> compute(SGNDArray<float64_t> C,
						   SGMatrix<float64_t> V0 = SGMatrix<float64_t>(NULL,0,0,false),
						   double eps=CMath::MACHINE_EPSILON,
						   int itermax=200) = 0;

		/** return the matrix V that best diagonalizes C */
		SGMatrix<float64_t> get_V()
		{
			return m_V;
		}

	protected:
		/** the matrix V that best diagonalizes C */
		SGMatrix<float64_t> m_V;

};

}


#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdio.h>
#include <assert.h>
#include <map>
#include <vector>
#include <string>

namespace shogun
{
class MethodProfiler
{
public :
    FILE* detailed_output = 0;
    int verbose;
    SGNDArray<float64_t> C;
    double t_start;
    double t_last_iteration;
    std::map<std::string, double> extra_output;
    int num_extra_params;
    
    MethodProfiler(const SGNDArray<float64_t>& C_, int verbose_, const char* file_name, std::vector<std::string> extra_output_cols = {}):
        verbose(verbose_), C(C_.clone()), t_start(clock()),
        t_last_iteration(t_start)
    {
        for (std::string s: extra_output_cols)
            extra_output[s] = 0;
        num_extra_params = extra_output.size();
        
        if (verbose == 2)
        {
            detailed_output = fopen(file_name, "w");
	    if (!detailed_output)
	    {
		    assert(!"Failed to open output file for method details");
	    }
            fprintf(detailed_output, "iteration,time,eps,quality_avg");
            for (int i = 0; i < C.dims[2]; i++)
                fprintf(detailed_output, ",quality_%d", i + 1);
            for (auto it: extra_output)
                fprintf(detailed_output, ",%s", it.first.c_str());
            fprintf(detailed_output, "\n");
            fflush(detailed_output);
        }
    }
    
    ~MethodProfiler()
    {
        if (detailed_output)
            fclose(detailed_output);
    }
    
    void iteration(int iter, const SGMatrix<float64_t>& V, double eps)
    {
        double t_last = t_last_iteration;
        double t_now = clock();
    
        if ((verbose == 1) || (verbose == 2))
        {
			printf("iteration %d done in %.3f s, eps %g", iter, (t_now - t_last) * 1.0 / CLOCKS_PER_SEC, eps);
            for (auto it: extra_output)
                printf(", %s %g", it.first.c_str(), it.second);
            printf("\n");
        }
	
        if (verbose == 2)
        {
            SGVector<float64_t> quality(C.dims[2]);
        
            for (int i = 0; i < C.dims[2]; i++)
            {
                SGMatrix<float64_t> A_diag =
                    SGMatrix<float64_t>::matrix_multiply(V,
                        SGMatrix<float64_t>::matrix_multiply(SGMatrix<float64_t>(C.get_matrix(i), V.num_cols, V.num_cols, false), V, false, true), false, false);
                
                double quality_off = 0;
                double quality_in = 0;
                
                for (int j = 0; j < A_diag.num_rows; j++)
                {
                    for (int k = 0; k < A_diag.num_cols; k++)
                    {
                        if (k == j) quality_in += std::abs(A_diag(j, k));
                        else quality_off += std::abs(A_diag(j, k));
                    }
                }
                
                quality[i] = quality_in / (quality_in + quality_off);
            }
            
            printf("iteration %d quality: ", iter);
            for (int i = 0; i < C.dims[2]; i++)
            {
                printf("%.3f ", quality[i]);
            }
            printf("\n");
            
            // iteration, time, eps, quality average, detailed quality list, detailed eps
            fprintf(detailed_output, "%d,%g,%g,%g",
                iter, (t_now - t_start) * 1.0 / CLOCKS_PER_SEC, eps,
                std::accumulate(quality.data(), quality.data() + C.dims[2], 0.0) / (double)C.dims[2]);
            
            for (int i = 0; i < C.dims[2]; i++)
                fprintf(detailed_output, ",%g", quality[i]);
            
            assert((int)extra_output.size() == num_extra_params);
            for (auto it: extra_output)
                fprintf(detailed_output, ",%g", it.second);
            fprintf(detailed_output, "\n");
            fflush(detailed_output);
        }
    }
};
}


#endif //APPROXJOINTDIAGONALIZER_H_
