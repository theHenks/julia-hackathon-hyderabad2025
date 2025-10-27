# Defining new Distributions - How to contribute to a julia package?

In our hackathon session, we want to develop a new distribution implementation in [DistributionsHEP.jl](https://github.com/JuliaHEP/DistributionsHEP.jl), namely the relativisitc Breit-Wigner distribution. To do so, we can use the exisitng code we developed during the tutorial part to test our implementation against our generic function implementation.
During the hackathon, we should guide our effort along the following list:

1. Get familiar with the `Distributions.jl` and `DistributionsHEP.jl` packages. Have a look into the `src/` folder, get familiar with the general package structure and how it is organized. 
2. Create your own fork of the repository and clone it to your computing setup. Add the package to your current environment using [Revise.jl](https://github.com/timholy/Revise.jl) with 
``` julia
] dev path/to/your/cloned/repo
```
3. Implement the relativistic Breit-Wigner distribution as a new type in `src/RelativisticBreitWigner.jl`. Implement all necessary methods to make it work with the `Distributions.jl` ecosystem. You can use any existing distribution implementation as a template. Make sure to implement at least the following methods:
   - `pdf`
   - `logpdf`
   - `rand`
   - `mean`
   - `var`
4. Write tests for your new distribution in the `test/` folder. You can use the existing tests as a template. Make sure to test all implemented methods and check that they behave as expected. In addition, check with an exisiting implementation of the relativistic Breit-Wigner distribution (e.g. your own generic function implementation from the tutorial part) that your implementation gives consistent results.
5. Document your new distribution in the `docs/` folder. Add a new section to the documentation that describes the relativistic Breit-Wigner distribution, its parameters, and how to use it.
6. Once you are satisfied with your implementation, tests, and documentation, create a pull request to the main `DistributionsHEP.jl` repository. Make sure to follow the contribution guidelines provided in the repository.
7. Celebrate your contribution to the Julia HEP ecosystem!


## Advanced task - Implementing convolution of distributions
As an advanced task, you can try to implement the convolution of two distributions in `DistributionsHEP.jl`. This is a more complex task that requires a deeper understanding of the `Distributions.jl` ecosystem and numerical methods. You can use the Breit-Wigner convoluted with a Gaussian as a test case for your implementation. Make sure to test your implementation thoroughly and document it properly. The method you need to dispatch on can be found [here](https://juliastats.org/Distributions.jl/stable/convolution/). For our fitting example and scenario, you can implement the convolution of the `RelativisticBreitWigner` distribution with the `Normal` distribution from `Distributions.jl`. Test the output of your convolution implementation against your own numerical convolution implementation. Try to implement the fit of the Z boson invariant mass spectrum using your new convolution distribution implementation.
Compare the results with the previous fits using your goodness-of-fit metrics.