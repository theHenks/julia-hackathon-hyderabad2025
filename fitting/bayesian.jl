using UnROOT
using FHist, StatsBase, Distributions, IntervalSets, LinearAlgebra
using Optimization, OptimizationLBFGSB, OptimizationOptimJL, QuadGK, ForwardDiff
using CairoMakie, LaTeXStrings
using Unitful, Measurements
using BAT, ValueShapes, DensityInterface

tfile = ROOTFile("data/zmass.root")

h = UnROOT.parseTH(tfile["Zmass"]; raw=false)
bin_width = only(unique(diff(h.binedges...)))

gaussian(x::AbstractFloat, p::NamedTuple{(:N, :μ, :σ)}) = p.N*pdf(Normal(p.μ, p.σ), x)
function relativistic_breit_wigner(m::AbstractFloat, p::NamedTuple{(:N, :M, :Γ)})
    let N = p.N, M = p.M, Γ = p.Γ, γ = √(M^2 * (M^2 + Γ^2)), k = (2*√(2)*M*Γ*γ) / (π * √(M^2 + γ))
        N * k / ((m^2 - M^2)^2 + M^2 * Γ^2)
    end
end

begin
    fig = Figure(size=(600, 400))
    ax = Axis(fig[1, 1]; 
                limits=((first(h.binedges[1]), last(h.binedges[1])), (-100, maximum(h.bincounts)*1.2)),
                xlabel=L"m_{Z} \; (\mathrm{GeV}/c^{2})", 
                ylabel=latexstring("Events / $bin_width \$\\mathrm{GeV}/c^{2}\$"),
                ylabelsize=18,
                xlabelsize=18,
                title=latexstring("\$Z_0\$ boson mass peak"),
                titlesize=20)
    errorbars!(ax, h, color=:black, whiskerwidth=6, label = "Candidate events")
    axislegend(ax; position = :rt)
    fig
end


function negative_hist_loglike(f_fit::Base.Callable, h::Hist1D{T}; fit_range::Tuple{T, T} = extrema(first(h.binedges))) where T<:AbstractFloat
    bin_edges = first(h.binedges)
    bin_centers = midpoints(bin_edges)
    fit_range_idxs = findall(in(Interval(fit_range...)), bin_centers)
    counts = h.bincounts[fit_range_idxs]
    bin_widths = diff(bin_edges)[fit_range_idxs]
    filter!(in(Interval(fit_range...)), bin_centers)
    # bin_ll(n, μ) =  μ - n * log(μ)
    # sum(bin_ll.(counts, f_fit.(bin_centers) .* bin_widths))
    bin_ll(x, bw, k) = logpdf(Poisson(bw * abs(f_fit(x)) + eps(T)), k)
    -sum(Base.Broadcast.broadcasted(bin_ll, bin_centers, bin_widths, counts))
end

# create loglikehood function: f_loglike(v) that can be evaluated for any set of v (fit parameter)
function f_loglike(f::Base.Callable; fit_range::Tuple{<:AbstractFloat, <:AbstractFloat} = extrema(first(h.binedges)))
    let h = h, fit_range = fit_range
        logfuncdensity(function (v)
            negative_hist_loglike(x -> f(x, v), h; fit_range=fit_range)
        end)
    end
end


# fit range and start values
fit_range = (30.0, 150.0)
v_prior_gaussian = distprod(N = Normal(sum(h.bincounts), sqrt(sum(h.bincounts))), μ = Normal(90.0, 5.0), σ = Weibull(10.0, 30.0))

loglike_init_gaussian = logdensityof(f_loglike(gaussian; fit_range=fit_range), mean(v_prior_gaussian))

# MCMC gaussian
posterior = PosteriorMeasure(f_loglike(gaussian), v_prior_gaussian)

ENV["JULIA_DEBUG"] = "Main"
samples = bat_sample(posterior, TransformedMCMC(proposal = RandomWalk(), nsteps = 10^5, nchains = 8, convergence=BrooksGelmanConvergence(threshold=10.0)))

# fit range and start values for Breit-Wigner
v_prior_bw = distprod(N = Normal(sum(h.bincounts), sqrt(sum(h.bincounts))), M = Normal(90.0, 5.0), Γ = Weibull(2.0, 5.0))

loglike_init_bw = logdensityof(f_loglike(relativistic_breit_wigner; fit_range=fit_range), mean(v_prior_bw))

# MCMC Breit-Wigner
posterior_bw = PosteriorMeasure(f_loglike(relativistic_breit_wigner), v_prior_bw)
samples_bw = bat_sample(posterior_bw, TransformedMCMC(proposal = RandomWalk(), nsteps = 10^5, nchains = 8))

let x_fit = first(fit_range):0.01:last(fit_range)
    fig = Figure(size=(600, 400))
    ax = Axis(fig[1, 1]; 
                limits=((first(h.binedges[1]), last(h.binedges[1])), (-100, maximum(h.bincounts)*1.3)),
                xlabel=L"m_{Z} \; (\mathrm{GeV}/c^{2})", 
                ylabel=latexstring("Events / $bin_width \$\\mathrm{GeV}/c^{2}\$"),
                ylabelsize=18,
                xlabelsize=18,
                title=latexstring("\$Z_0\$ boson mass peak"),
                titlesize=20)
    errorbars!(ax, h, color=:black, whiskerwidth=6, label = "Candidate events")
    band!(ax, x_fit, x -> gaussian(x, v_best_gaussian_errx10...)*bin_width; color=:orange, alpha=0.2, label="Gaussian fit (1σ Err x10)")
    lines!(ax, x_fit, x -> gaussian(x, v_best_gaussian_errx10...)*bin_width; color=:orange, linewidth=2, label="Gaussian fit (1σ Err x10)")
    band!(ax, x_fit, x -> relativistic_breit_wigner(x, v_best_bw_errx10...)*bin_width; color=:darkblue, alpha=0.2, label="Breit-Wigner fit (1σ Err x10)")
    lines!(ax, x_fit, x -> relativistic_breit_wigner(x, v_best_bw_errx10...)*bin_width; color=:darkblue, linewidth=2, label="Breit-Wigner fit (1σ Err x10)")
    axislegend(ax; position = :rt, merge=true, unique=true)
    fig
end