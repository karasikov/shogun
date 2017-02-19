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

#ifndef MMD_H_
#define MMD_H_

#include <utility>
#include <memory>
#include <functional>
#include <shogun/statistical_testing/TwoSampleTest.h>

namespace shogun
{

class CKernel;
class CCustomDistance;
template <typename> class SGVector;
template <typename> class SGMatrix;

namespace internal
{

class MaxTestPower;
class MaxXValidation;
class WeightedMaxTestPower;

}

enum EStatisticType
{
	ST_UNBIASED_FULL,
	ST_UNBIASED_INCOMPLETE,
	ST_BIASED_FULL
};

enum EVarianceEstimationMethod
{
	VEM_DIRECT,
	VEM_PERMUTATION
};

enum ENullApproximationMethod
{
	NAM_PERMUTATION,
	NAM_MMD1_GAUSSIAN,
	NAM_MMD2_SPECTRUM,
	NAM_MMD2_GAMMA
};

enum EKernelSelectionMethod
{
	KSM_MEDIAN_HEURISTIC,
	KSM_MAXIMIZE_MMD,
	KSM_MAXIMIZE_POWER,
	KSM_MAXIMIZE_XVALIDATION,
	KSM_AUTO
};

class CMMD : public CTwoSampleTest
{
	using operation=std::function<float32_t(SGMatrix<float32_t>)>;
	friend class internal::MaxTestPower;
	friend class internal::WeightedMaxTestPower;
	friend class internal::MaxXValidation;
public:
	CMMD();
	virtual ~CMMD();

	void add_kernel(CKernel *kernel);
	void select_kernel(EKernelSelectionMethod kmethod=KSM_AUTO,
		bool weighted_kernel=false, float64_t train_test_ratio=1.0,
		index_t num_run=10, float64_t alpha=0.05);

	virtual float64_t compute_statistic() override;
	virtual float64_t compute_variance();

	virtual SGVector<float64_t> sample_null() override;

	void use_gpu(bool gpu);
	void cleanup();

	void set_statistic_type(EStatisticType stype);
	const EStatisticType get_statistic_type() const;

	void set_variance_estimation_method(EVarianceEstimationMethod vmethod);
	const EVarianceEstimationMethod get_variance_estimation_method() const;

	void set_num_null_samples(index_t null_samples);
	const index_t get_num_null_samples() const;

	void set_null_approximation_method(ENullApproximationMethod nmethod);
	const ENullApproximationMethod get_null_approximation_method() const;

	virtual const char* get_name() const;
protected:
	virtual const operation get_direct_estimation_method() const=0;
	virtual const float64_t normalize_statistic(float64_t statistic) const=0;
	virtual const float64_t normalize_variance(float64_t variance) const=0;
	bool use_gpu() const;
private:
	struct Self;
	std::unique_ptr<Self> self;
	virtual std::pair<float64_t, float64_t> compute_statistic_variance();
	virtual std::pair<SGVector<float64_t>, SGMatrix<float64_t>> compute_statistic_and_Q();
};

}
#endif // MMD_H_
