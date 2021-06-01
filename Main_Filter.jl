using LinearAlgebra
#using Pkg
#Pkg.add("Distributions")
using Distributions
#Pkg.add("Combinatorics")
import Combinatorics: combinations #combinations(a,n): returns all combinations of n elements of indexable object a;
#Pkg.add("JLD")
using JLD

include("configure_sched_file.jl")

function constraints(x::Vector; Nx = 42, Ny = 42, distance_limit = 20)
	c1 = [-l+1 for l in x]
	c2 = [l-Nx*Ny for l in x]
	c3 = [-distance(ls[1], ls[2])+distance_limit for ls in collect(combinations(x,2))]
	return [c1; c2; c3]
end

function get_l(i, j; Nx = 42, Ny = 42)
	return l = floor(Int,(j-1)*Nx + i)
end

function get_coordinate(l::Int; Nx = 42)
	j = ceil(Int, l/Nx)
	i = l - (j-1)*Nx
	return [i, j]
end

function distance(l1::Int,l2::Int)
	vec1 = get_coordinate(l1)
	vec2 = get_coordinate(l2)
	return norm(vec1-vec2)
end

function wells_inside_domain(well_dict::Dict, x, c)
	N = length(well_dict) #number of wells
	c_value = c(convert.(Int, round.(x)))
	if sum(c_value[1:2*N].>= 0) == 0 #then there are no boundary constraints violations
		return true		
	else #if an element is >0 then there's a well outside the reservoir boundary and Eclipse won't run
		return 	false
	end
end

function covariance_matrix_adaptation(c, k_max, σ)
	directory = pwd()

	#Read initial wells' coordinates from SCHED file that is placed in the same directory as the current code
	well_dict = read_sched_file("base.SCHED", directory)
	global wells_name =  [i for i in keys(well_dict)]

	#Temporary solution before the decoder is implemented:
	#Sample just the x coordinate of each well (2D reservoir)
	x = [get_l(i[2][1],i[2][2]) for i in well_dict]
	x_history = x #initialize vector x_history
    ys_history = [] #initialize variable

	m = 4 + floor(Int, 3*log(length(x)))
	m_elite = div(m,2)
	μ, n = copy(x), length(x) #n is the number of wells
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

	#Create m folders inside current directory
	for i = 1:m
		rm(string(i), force = true, recursive=true) #delete folder if it exists
		mkdir(string(i))
		mkpath(string(i) * "\\AA222_EclipseFile")
		file_directory = directory * "\\" * string(i)
		#Copy the folder "A222_EclipseFile" and the file "FEVAL.py" to each folder
		#Couldn't find a function that copies and pastes the entire folder.... So I'm doing it file by file
		cp(directory * "\\AA222_EclipseFile\\AA222.data", file_directory * "\\AA222_EclipseFile\\AA222.data", force = true)
		cp(directory * "\\AA222_EclipseFile\\base.sched", file_directory * "\\AA222_EclipseFile\\base.sched", force = true)
		cp(directory * "\\AA222_EclipseFile\\PERMX.IN", file_directory * "\\AA222_EclipseFile\\PERMX.IN", force = true)
		cp(directory * "\\FEVAL.py", file_directory  * "\\FEVAL.py", force = true)
	end

	NPV = zeros(m) #initialize vector

	for k in 1 : k_max
		println(k)
		
		P = MvNormal(μ, σ^2*Σ)
		xs = [rand(P) for i in 1 : m]

		#Configure each SHED file with new wells coordinates and place it in its respective folder
		#And call the Python script that runs ECLIPSE simulator inside each individual folder
		#And read the NPV.text file inside each folder and store the values in the NPV variable

		for i_individual = 1:m
			#println("xs[i_individual] = $(xs[i_individual])")
			#println("wells_inside_domain? $(wells_inside_domain(well_dict, xs[i_individual], c))")

			if wells_inside_domain(well_dict, xs[i_individual], c) == true
				#From x variable build well_dict:
				for i_well in keys(well_dict)
					i_well_index = findfirst((x -> x==i_well), wells_name)
					i,j = get_coordinate(convert(Int, round(xs[i_individual][i_well_index], digits=0)))
					well_dict[i_well][1:2] = [i,j] #change only the x and y positions on this test case
				end
				file_directory = directory * "\\" * string(i_individual) * "\\" * "AA222_EclipseFile"
				write_sched_file("base.sched", well_dict, file_directory)

				pythondirectory = directory * "\\" * string(i_individual)
				cd(pythondirectory)
				mycommand = `python FEVAL.py`
				run(mycommand)
				cd(directory)
				
				aux = open(f->read(f, String), directory * "\\" * string(i_individual) * "\\AA222_EclipseFile\\NPV.txt")
				NPV[i_individual] = parse(Float64, aux)
			else #then the boundary constraints are violated
				NPV[i_individual] = 0
			end
		end

        constraint = [sum(max.(c(convert.(Int, round.(x))),0).^2) for x in xs] #Quadratic penalty function
		ys = [(constraint[i],-NPV[i],i) for i in 1:m]
        #println("ys = $ys")
		is = sortperm(ys, by = x->(x[1],x[2])) # best to worst
        #println("is = $is")

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

		x_history = [x_history convert.(Int, round.(μ))]
        push!(ys_history, ys)

	end
	return x_history, ys_history
end

function optimize(c, k_max, σ)
	return (x_history, ys_history) = covariance_matrix_adaptation(c, k_max, σ)
end

σ_vetor = [500]#[1, 5, 10, 50, 100, 500]
for i_σ in 1:length(σ_vetor)
    σ = σ_vetor[i_σ]
    cd(dirname(@__FILE__)) #change location to current directory
    (x_history, ys_history) = optimize(constraints, 100, σ)
    save("resultsDIST20_filter_sigma$σ.jld", "x_history", x_history, "ys_history", ys_history)

end