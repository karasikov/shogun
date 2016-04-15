/*
 * Copyright (c) The Shogun Machine Learning Toolbox
 * Written (w) 2014 Wu Lin
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
 *
 * the reference paper is
 * Mohammad Emtiyaz Khan, Aleksandr Y. Aravkin, Michael P. Friedlander, Matthias Seeger
 * Fast Dual Variational Inference for Non-Conjugate Latent Gaussian Models. ICML2013
 *
 */

#ifndef _KLDUALINFERENCEMETHOD_H_
#define _KLDUALINFERENCEMETHOD_H_

#include <shogun/lib/config.h>
#include <shogun/optimization/lbfgs/lbfgs.h>
#include <shogun/machine/gp/KLInferenceMethod.h>
#include <shogun/machine/gp/DualVariationalGaussianLikelihood.h>

namespace shogun
{

/** @brief The dual KL approximation inference method class.
 *
 * This inference process is described in the reference paper
 * Mohammad Emtiyaz Khan, Aleksandr Y. Aravkin, Michael P. Friedlander, Matthias Seeger
 * Fast Dual Variational Inference for Non-Conjugate Latent Gaussian Models. ICML2013
 *
 * The idea is to optimize the log marginal likelihood
 * with equality constaints (primal problem) by solving the Lagrangian dual problem.
 * The equality constaints are:
 * \f[
 * h = \mu, \rho = \sigma^2 = diag(\Sigma)
 * \f], where h and \f$\rho\f$ are auxiliary variables, \f$\mu\f$ and \f$\sigma^2\f$ are variational variables,
 * and \f$\Sigma\f$ is an approximated posterior covariance matrix.
 * The equality constaints are variational mean constaint (\f$\mu\f$) and variational variance constaint (\f$\sigma^2\f$).
 *
 * For detailed information, please refer to the paper.
 */
class CKLDualInferenceMethod: public CKLInferenceMethod
{
public:
	/** default constructor */
	CKLDualInferenceMethod();

	/** constructor
	 *
	 * @param kernel covariance function
	 * @param features features to use in inference
	 * @param mean mean function
	 * @param labels labels of the features
	 * @param model Likelihood model to use
	 */
	CKLDualInferenceMethod(CKernel* kernel, CFeatures* features,
			CMeanFunction* mean, CLabels* labels, CLikelihoodModel* model);

	virtual ~CKLDualInferenceMethod();

	/** returns the name of the inference method
	 *
	 * @return name KLDualInferenceMethod
	 */
	virtual const char* get_name() const { return "KLDualInferenceMethod"; }

	/** return what type of inference we are
	 *
	 * @return inference type KL_DUAL
	 */
	virtual EInferenceType get_inference_type() const { return INF_KL_DUAL; }

	/** helper method used to specialize a base class instance
	 *
	 * @param inference inference method
	 * @return casted CKLDualInferenceMethod object
	 */
	static CKLDualInferenceMethod * obtain_from_generic(CInferenceMethod* inference);

	/** get alpha vector
	 *
	 * @return vector to compute posterior mean of Gaussian Process:
	 *
	 * \f[
	 * \mu = K\alpha+mean_f
	 * \f]
	 *
	 * where \f$\mu\f$ is the mean and \f$K\f$ is the prior covariance matrix.
	 */
	virtual SGVector<float64_t> get_alpha();

	/** get diagonal vector
	 *
	 * @return diagonal of matrix used to calculate posterior covariance matrix:
	 *
	 * \f[
	 * Cov = (K^{-1}+sW^{2})^{-1}
	 * \f]
	 *
	 * where \f$Cov\f$ is the posterior covariance matrix, \f$K\f$ is the prior
	 * covariance matrix, and \f$sW\f$ is the diagonal vector.
	 */
	virtual SGVector<float64_t> get_diagonal_vector();

	/** set variational likelihood model
	 *
	 * @param mod model to set
	 */
	void set_model(CLikelihoodModel* mod);


	/** set L-BFGS parameters
	 * For details please see shogun/optimization/lbfgs/lbfgs.h
	 * @param m The number of corrections to approximate the inverse hessian matrix.
	 * Default value is 100.
	 * @param max_linesearch The maximum number of trials to do line search for each L-BFGS update.
	 * Default value is 1000.
	 * @param linesearch The line search algorithm.
	 * Default value is using the backtracking with the strong Wolfe condition line search
	 * @param max_iterations The maximum number of iterations for L-BFGS update.
	 * Default value is 1000.
	 * @param delta Delta for convergence test based on the change of function value.
	 * Default value is 0.
	 * @param past Distance for delta-based convergence test.
	 * Default value is 0.
	 * @param epsilon Epsilon for convergence test based on the change of gradient.
	 * Default value is 1e-5
	 * @param min_step The minimum step of the line search.
	 * The default value is 1e-20
	 * @param max_step The maximum step of the line search.
	 * The default value is 1e+20
	 * @param ftol A parameter used in Armijo condition.
	 * Default value is 1e-4
	 * @param wolfe A parameter used in curvature condition.
	 * Default value is 0.9
	 * @param gtol A parameter used in Morethuente linesearch to control the accuracy.
	 * Default value is 0.9
	 * @param xtol The machine precision for floating-point values.
	 * Default value is 1e-16.
	 * @param orthantwise_c Coeefficient for the L1 norm of variables.
	 * This parameter should be set to zero for standard minimization problems.
	 * Setting this parameter to a positive value activates
	 * Orthant-Wise Limited-memory Quasi-Newton (OWL-QN) method. Default value is 0.
	 * @param orthantwise_start Start index for computing L1 norm of the variables.
	 * This parameter is valid only for OWL-QN method. Default value is 0.
	 * @param orthantwise_end End index for computing L1 norm of the variables.
	 * Default value is 1.
	 */
	virtual void set_lbfgs_parameters(int m = 100,
		int max_linesearch = 1000,
		int linesearch = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE,
		int max_iterations = 1000,
		float64_t delta = 0.0,
		int past = 0,
		float64_t epsilon = 1e-5,
		float64_t min_step = 1e-20,
		float64_t max_step = 1e+20,
		float64_t ftol = 1e-4,
		float64_t wolfe = 0.9,
		float64_t gtol = 0.9,
		float64_t xtol = 1e-16,
		float64_t orthantwise_c = 0.0,
		int orthantwise_start = 0,
		int orthantwise_end = 1);

protected:

	/** The number of corrections to approximate the inverse hessian matrix.*/
	int m_m;

	/** The maximum number of trials to do line search for each L-BFGS update.*/
	int m_max_linesearch;

	/** The line search algorithm.*/
	int m_linesearch;

	/** The maximum number of iterations for L-BFGS update.*/
	int m_max_iterations;

	/** Delta for convergence test based on the change of function value.*/
	float64_t m_delta;

	/** Distance for delta-based convergence test.*/
	int m_past;

	/** Epsilon for convergence test based on the change of gradient.*/
	float64_t m_epsilon;

	/** The minimum step of the line search.*/
	float64_t m_min_step;

	/** The maximum step of the line search.*/
	float64_t m_max_step;

	/** A parameter used in Armijo condition.*/
	float64_t m_ftol;

	/** A parameter used in curvature condition.*/
	float64_t m_wolfe;

	/** A parameter used in Morethuente linesearch to control the accuracy.*/
	float64_t m_gtol;

	/** The machine precision for floating-point values.*/
	float64_t m_xtol;

	/** Coeefficient for the L1 norm of variables.*/
	float64_t m_orthantwise_c;

	/** Start index for computing L1 norm of the variables.*/
	int m_orthantwise_start;

	/** End index for computing L1 norm of the variables.*/
	int m_orthantwise_end;

	/** compute the gradient wrt variational parameters
	 * given the current variational parameters (mu and s2)
	 *
	 * @return gradient of negative log marginal likelihood
	 */
	virtual void get_gradient_of_nlml_wrt_parameters(SGVector<float64_t> gradient){};

	/** this method is used to dynamic-cast the likelihood model, m_model,
	 * to dual variational likelihood model.
	 */
	virtual CDualVariationalGaussianLikelihood* get_dual_variational_likelihood() const;

	/** check the provided likelihood model supports dual variational inference or not
	 * @param mod the provided likelihood model
	 *
	 * @return whether the provided likelihood model supports dual variational inference or not
	 */
	virtual void check_dual_inference(CLikelihoodModel* mod) const;

	/** update covariance matrix of the approximation to the posterior */
	virtual void update_approx_cov();

	/** update alpha matrix */
	virtual void update_alpha();

	/** update cholesky matrix */
	virtual void update_chol();

	/** update matrices which are required to compute negative log marginal
	 * likelihood derivatives wrt hyperparameter
	 */
	virtual void update_deriv();

	/** the helper function to compute
	 * the negative log marginal likelihood
	 *
	 * @return negative log marginal likelihood
	 */
	virtual float64_t get_negative_log_marginal_likelihood_helper();

	/** pre-compute the information for optimization.
	 * This function needs to be called before calling
	 * get_negative_log_marginal_likelihood_wrt_parameters()
	 * and/or
	 * get_gradient_of_nlml_wrt_parameters(SGVector<float64_t> gradient)
	 *
	 * @return true if precomputed parameters are valid
	 */
	virtual bool precompute();

	/** compute matrices which are required to compute negative log marginal
	 * likelihood derivatives wrt  hyperparameter in cov function
	 * Note that
	 * get_derivative_wrt_inference_method(const TParameter* param)
	 * and
	 * get_derivative_wrt_kernel(const TParameter* param)
	 * will call this function
	 *
	 * @param dK the gradient related to cov
	 *
	 * @return the gradient wrt hyperparameter related to cov
	 */
	virtual float64_t get_derivative_related_cov(SGMatrix<float64_t> dK);

	/** Using L-BFGS to estimate posterior parameters */
	virtual float64_t optimization();

	/** compute the objective value for LBFGS optimizer
	 *
	 * The mathematical equation (equation 24 in the paper) is defined as below
	 * \f[
	 * min_{\lambda \in S}{0.5*[(\lambda-y)^TK(\lambda-y)-log(det(A_{\lambda}))]-mean_{f}^T(\lambda-y)+\sum_{i=1}^{n}{Fenchel_i{(\lambda)}}}
	 * \f]
	 * where S is the feasible set defined for \f$\lambda\f$, K comes from covariance function,
	 * \f$mean_f\f$ comes from mean function, \f$\lambda\f$ is the dual parameter,
	 * y are data labels, n is the number point,
	 * \f$A_{\lambda}=K^{-1}+diag(\lambda)\f$,
	 * and \f$Fenchel_i{(\lambda)}=Fenchel_i{(\alpha,\lambda)}\f$ since \f$\alpha\f$ is implicitly defined by \f$\lambda\f$
	 *
	 * Note that S and \f$Fenchel_i{(\lambda)}\f$ are specified by the data modeling distribution,
	 * which are implemented in dual variational likelihood class.
	 *
	 */
	virtual float64_t get_dual_objective_wrt_parameters();

	/** compute the gradient of the objective function for LBFGS optimizer
	  * The mathematical equation (equation 25 in the paper) is defined as below
	  * \f[
	  * 0.5*[2*K(\lambda-y)-diag(A_{\lambda}^{-1})]-mean_{f}+\sum_{i=1}^{n}{\nabla Fenchel_i{(\lambda)}}
	  * \f]
	  * where \f$A_{\lambda}=K^{-1}+diag(\lambda)\f$,
	  * K comes from covariance function,
	  * \f$mean_f\f$ comes from mean function,
	  * \f$\lambda\f$ is the dual parameter,
	  * y are data labels, n is the number point,
	  * and \f$\nabla Fenchel_i{(\lambda)}\f$ is the gradient of \f$Fenchel_i{(\lambda)}\f$ wrt to \f$\lambda\f$
	  *
	  *Note that \f$\nabla Fenchel_i{(\lambda)}\f$ are specified by the data modeling distribution,
	  *which are implemented in dual variational likelihood class.
	  */
	virtual void get_gradient_of_dual_objective_wrt_parameters(SGVector<float64_t> gradient);

private:
	void init();

	/** square root of noise matrix W */
	SGVector<float64_t> m_sW;

	/** noise Matrix
	 * Note that it is also the dual parameter, lambda
	 */
	SGVector<float64_t> m_W;

	/** the gradient of the variational expection wrt sigma2
	 * sigma2=diag(m_Sigma)
	 */
	SGVector<float64_t> m_dv;

	/** the gradient of the variational expection wrt mu*/
	SGVector<float64_t> m_df;

	/** the Matrix V, where
	 * L'*V=diag(sW)*K
	 * Note that L' is a lower triangular matrix
	 */
	SGMatrix<float64_t> m_V;
	/**
	 * whether the lambda (m_W) is valid or not.
	 * In other word, whether it is feasible or not
	 * according to the dual variational likelihood
	 */
	bool m_is_dual_valid;

	/** compute the negative log marginal likelihood
	 * given the current variational parameters (mu and s2)
	 *
	 * @param alpha
	 * @param mu
	 * @param L
	 *
	 * @return negative log marginal likelihood
	 */

	float64_t get_nlml_wrapper(SGVector<float64_t> alpha, SGVector<float64_t> mu, SGMatrix<float64_t> L);

	/** helper function is passed to the LBFGS API
	 * to compute objective value and gradient
	 * Note that this function should be static
	 * */
	static float64_t evaluate(void *obj,
		const float64_t *parameters,
		float64_t *gradient, const int dim,
		const float64_t step);

	/** helper function is passed to the LBFGS API
	 * to adjust step size based on the feasible set S
	 * defined in dual variational likelihood.
	 * Note that this function should be static
	 * */
	static float64_t adjust_step(void *obj,
		const float64_t *parameters,
		const float64_t *direction,
		const int dim, const float64_t step);

};
}
#endif /* _KLDUALINFERENCEMETHOD_H_ */
