using LinearAlgebra
using Distributions

include("helpers.jl")
include("simple.jl")

#include("cmaes.jl")

function penalty_funtion(choice)
	fun = if choice == "quadratic"
		function pquadratic(x, c)
			return sum(max.(c(x),0).^2)
		end
	elseif choice == "count"
		function pcount(x, c)
			return sum(sum([c(x) .> 0])) #gives the total number of constraint violations
		end
	end
	return fun
end

function ∇penalty_funtion(choice)
	fun = if choice == "quadratic"
		function ∇pquadratic(x, c)
			#TO BE IMPLEMENTED
			println("TO BE IMPLEMENTED")
			return 0
		end
	elseif choice == "count"
		function ∇pcount(x, c)
			return 0
		end
	end
	return fun
end

function covariance_matrix_adaptation(f, x, k_max;σ = 1.0)
	x_history = x

	m = 4 + floor(Int, 3*log(length(x)))
	m_elite = div(m,2)
	μ, n = copy(x), length(x)
	ws = log((m+1)/2) .- log.(1:m)
	ws[1:m_elite] ./= sum(ws[1:m_elite])
	μ_eff = 1 / sum(ws[1:m_elite].^2)
	cσ = (μ_eff + 2)/(n + μ_eff + 5)
	dσ = 1 + 2max(0, sqrt((μ_eff-1)/(n+1))-1) + cσ
	cΣ = (4 + μ_eff/n)/(n + 4 + 2μ_eff/n)
	c1 = 2/((n+1.3)^2 + μ_eff)
	cμ = min(1-c1, 2*(μ_eff-2+1/μ_eff)/((n+2)^2 + μ_eff))
	ws[m_elite+1:end] .*= -(1 + c1/cμ)/sum(ws[m_elite+1:end])
	E = n^0.5*(1-1/(4n)+1/(21*n^2))
	pσ, pΣ, Σ = zeros(n), zeros(n), Matrix(1.0I, n, n)
	for k in 1 : k_max
		P = MvNormal(μ, σ^2*Σ)
		xs = [rand(P) for i in 1 : m]
		ys = [f(x) for x in xs]
		is = sortperm(ys) # best to worst
		# selection and mean update
		δs = [(x - μ)/σ for x in xs]
		δw = sum(ws[i]*δs[is[i]] for i in 1 : m_elite)
		μ += σ*δw
		# step-size control
		C = Σ^-0.5
		pσ = (1-cσ)*pσ + sqrt(cσ*(2-cσ)*μ_eff)*C*δw
		σ *= exp(cσ/dσ * (norm(pσ)/E - 1))
		# covariance adaptation
		hσ = Int(norm(pσ)/sqrt(1-(1-cσ)^(2k)) < (1.4+2/(n+1))*E)
		pΣ = (1-cΣ)*pΣ + hσ*sqrt(cΣ*(2-cΣ)*μ_eff)*δw
		w0 = [ws[i]≥0 ? ws[i] : n*ws[i]/norm(C*δs[is[i]])^2 for i in 1:m]
		Σ = (1-c1-cμ) * Σ +
		c1*(pΣ*pΣ' + (1-hσ) * cΣ*(2-cΣ) * Σ) +
		cμ*sum(w0[i]*δs[is[i]]*δs[is[i]]' for i in 1 : m)
		Σ = triu(Σ)+triu(Σ,1)' # enforce symmetry

		x_history = [x_history μ]
	end
	return x_history
end

function optimize(f, ∇f, x0, n_max, prob)
	return x_history = covariance_matrix_adaptation(f, x0, n_max;σ = 1.0)
end

x_history = optimize(rosenbrock, rosenbrock_gradient, rosenbrock_init(), 20, "simple1")
