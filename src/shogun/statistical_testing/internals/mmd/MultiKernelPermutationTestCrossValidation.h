/*
 * Copyright (c) The Shogun Machine Learning Toolbox
 * Written (w) 2016 Soumyajit De
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies,
 * either expressed or implied, of the Shogun Development Team.
 */

#ifndef MULTIKERNEL_PERMUTATION_TEST_CROSS_VALIDATION
#define MULTIKERNEL_PERMUTATION_TEST_CROSS_VALIDATION

#include <memory>
#include <shogun/lib/common.h>
#include <shogun/lib/SGMatrix.h>

namespace shogun
{

enum class EStatisticType;
class CCustomDistance;

namespace internal
{

class KernelManager;

namespace mmd
{

/**
 * @brief class that runs cross-validation test for MMD for multiple kernels.
 */
class MultiKernelPermutationTestCrossValidation
{
public:
	MultiKernelPermutationTestCrossValidation(index_t nx, index_t ny, index_t null_samples, EStatisticType type);
	~MultiKernelPermutationTestCrossValidation();

	void operator()(const KernelManager& km);

	void set_num_runs(index_t nr);
	void set_num_folds(index_t nf);
	void set_alpha(float64_t alp);
	void set_measure_matrix(SGMatrix<float64_t> measures);
	void set_distance(CCustomDistance* distance);
private:
	struct terms_t;

	void add_term(terms_t&, float64_t kernel, index_t i, index_t j);
	float64_t compute_mmd(terms_t&);

	const index_t n_x;
	const index_t n_y;
	const index_t num_null_samples;
	const EStatisticType stype;
	index_t num_runs;
	index_t num_folds;
	float64_t alpha;
	SGMatrix<float64_t> rejections;
	std::shared_ptr<CCustomDistance> m_distance;
};

}

}

}

#endif // MULTIKERNEL_PERMUTATION_TEST_CROSS_VALIDATION
