# Likelihood and fitting

In this tutorial, you will explore fitting continuous distributions to a data sample and parameter estimation.

You will be studying a sample that resembles the Z boson mass peak using [CMS data](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookHowToFit). The file contains a histogram of the invariant mass of reconstructed Z boson candidates. You can find it in `data/zmass.root` in the `gihtub` folder. 

## Details about the fits
To fit a generator-level Z peak a Breit-Wigner fit makes sense. However, reconstructed-level Z peaks have many detector resolutions that smear the Z mass peak. If the detector resolution is relatively poor, then it is usually good enough to fit a gaussian (since the gaussian detector resolution will overwhelm the inherent Breit-Wigner shape of the peak). If the detector resolution is fairly good, then another option is to fit a Breit-Wigner (for the inherent shape) convoluted with a gaussian (to describe the detector effects).This is in the "no-background" case. If you have backgrounds in your sample (Drell-Yan, cosmics, etc...), and you want to do the fit over a large mass range, then another function needs to be included to take care of this - an exponential is commonly used (cf. adavanced tutorial).
In this tutorial, we want to fit the data with 3 different models using our own likehood implementations:
1. [Gaussian](https://en.wikipedia.org/wiki/Normal_distribution):
$
    G(x ; \mu, \sigma)=\frac{1}{\sqrt{2 \pi} \sigma} \exp \left[-\frac{(x-\mu)^2}{2 \sigma^2}\right]
$
2. [Relativistic Breit-Wigner](https://en.wikipedia.org/wiki/Relativistic_Breit%E2%80%93Wigner_distribution):
$
    B(m; M, \Gamma)=\frac{k}{\left(E^2-M^2\right)^2+M^2 \Gamma^2}
$

where $k$ is defined as
$
    k=\frac{2 \sqrt{2} M \Gamma \gamma}{\pi \sqrt{M^2+\gamma}}, \quad \gamma=\sqrt{M^2\left(M^2+\Gamma^2\right)}
$

3. [Breit-Wigner convoluted with a Gaussian](https://en.wikipedia.org/wiki/Voigt_profile):
$
P(m)=\int B\left(m^{\prime} ; M, \Gamma\right) \cdot G\left(m-m^{\prime} ; \mu, \sigma\right) d m^{\prime}
$

## Exercise 1 - Loading the data

Open the provided ROOT file and extract the histogram of the invariant mass of reconstructed Z boson candidates. To open a ROOT file in Julia, you can use the [UnROOT.jl](https://github.com/JuliaHEP/UnROOT.jl).
The resulting histogram is parsed as a [`FHist`](https://github.com/Moelf/FHist.jl). Plot the histogram using [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl)


## Exercise 2 - MLE fitting

Fit the histogram of the invariant mass of reconstructed Z boson candidates with two of the three models mentioned above: Gaussian and Relativistic Breit-Wigner. Briefly explain your choice of the initial parameters for each fit. Start using a frequentist approach (e.g. maximum likelihood estimation) and then try a Bayesian approach (e.g. MCMC sampling). Compare the results of both approaches.

### Frequentist approach
- Implement the negative log-likelihood functions which can be evaluated given the histogram and a generic model function. 
- Define and implement the model functions for the three models mentioned above and set starting values for the fit parameters.
- Use an optimization package (e.g. [Optimization.jl](https://docs.sciml.ai/Optimization/stable/optimization_packages/optimization/) ) to minimize the negative log-likelihood functions and extract the best-fit parameters. Feel free to play with different optimizers and settings, e.g. LBFGS-B or Minuit. Try to activate the automatic differentiation if supported by the optimizer. 
- Extract the uncertainties on the fit parameters from the covariance matrix. To do so, you can compute the Hessian matrix of the negative log-likelihood function at the best-fit point and then invert it to obtain the covariance matrix using an automatic differentiation package (e.g. [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) ). 
- Plot the histogram of the invariant mass of reconstructed Z boson candidates along with the best-fit models obtained from the MLE fits. Include a $1 \sigma$ error band around the best-fit models using the uncertainties on the fit parameters.

### Bayesian approach
- Implement the log-likelihood functions which can be evaluated given the histogram and a generic model function. You can reuse the model and likehood functions defined for the frequentist approach.
- Make yourself familiar with the [BAT.jl](https://bat.github.io/BAT.jl/stable/) package for Bayesian analysis in Julia. Try to understand how to define prior distributions, posterior measures, and how to perform MCMC sampling. Choose a set of prior distributions for the fit parameters and justify your choices. Try out different prior distributions and see how they affect the results.
- Use `BAT.jl` to perform MCMC sampling of the posterior distributions of the fit parameters given the histogram data. Make sure to run the MCMC chains for a sufficient number of steps to ensure convergence and good mixing. In case no convergence is reached, try to tune the MCMC settings (e.g. step size, number of chains, etc...) or argue why convergence is not reached.
- Extract the best-fit parameters and their uncertainties from the posterior distributions. Extract the mean, median, mode and credible intervals (e.g. 68%, 95.5% and 99.7% credible interval) for each fit parameter. 
- Plot the histogram of the invariant mass of reconstructed Z boson candidates along with credibility intervals as well as the global mode and median. 

## Exercise 3 - Goodness of Fit

Evaluate the goodness of fit for each model and fitting approach. You can use any statistical test of your choice, e.g. 
- (Posterior predictiv) p-value
- Reduced $\Chi^2$
Compare the goodness of fit results for each model and fitting approach. Which model and approach provides the best fit to the data? Justify your answer based on the goodness of fit results. 