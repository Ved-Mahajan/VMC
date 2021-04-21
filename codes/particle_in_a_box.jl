### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ 7a0128c0-a13f-11eb-1098-ef9c00274a0d
using Plots

# ╔═╡ 7c5e40c4-f84e-4527-a87e-c6ac395f5d9e
function ansatz(x,α)
	return ifelse(abs(x) > 1,1 - abs(x)^α,0)
end;

# ╔═╡ 007bb513-b173-4723-9d8b-2e858ae54dc6
function probability_density(x,α)
	return (ansatz(x,α)^2)
end;

# ╔═╡ dcc6bdb1-daa8-44f1-a6ec-285ac1916179
function E_local(x,α)
	return ifelse(abs(x)<1,(0.5*α*(α - 1)*abs(x)^(α - 2))/(1 - abs(x)^α),0)
end;

# ╔═╡ 99eff10f-5b1e-42d9-9009-56e5aeb821c2
function metropolis(N,α)
	L = 2
	x = rand(-1:0.001:1)
	E = 0
	E2 = 0
	Eln_avg = 0
	ln_avg = 0
	rejection_ratio = 0
	#start the loop
	for i in 1:N
		x_t = x + 0.01*rand(-L:0.000001:L)
		if probability_density(x_t,α) >= probability_density(x,α)
			x = x_t
		else
			dummy = rand()
			if dummy < probability_density(x_t,α)/probability_density(x,α)
				x = x_t
			else
				rejection_ratio += 1/N
			end
		end
		
		E += E_local(x,α)/N
		E2 += E_local(x,α)^2/N
		Eln_avg += (E_local(x,α)*((-1*log(abs(x))*abs(x)^α)/(1 - abs(x)^α)))/(N)
		ln_avg += ((-1*log(abs(x))*abs(x)^α)/(1 - abs(x)^α))/N
	end
	return E,E2,Eln_avg,ln_avg,rejection_ratio
end;

# ╔═╡ cbfdf51f-721e-4172-b7d7-1911bf666163
function iterate_alpha()
	α = 1.1
	alpha_iter = 50
	N = 500
	random_walkers = 200
	γ = 0.01
	energy = Vector{Float64}(undef,alpha_iter)
	alpha = Vector{Float64}(undef,alpha_iter)
	variance = Vector{Float64}(undef,alpha_iter)
	for i in 1:alpha_iter
		E = 0
		E2 = 0
		dE_dalpha = 0
		Eln = 0
		ln = 0
		rejection_ratio = 0
		for j in 1:random_walkers
			E_met, E2_met, Eln_met, ln_met, _ = metropolis(N,α)
			E += E_met/random_walkers
			E2 += E2_met/random_walkers
			Eln += Eln_met/random_walkers
			ln += ln_met/random_walkers
			# rejection_ratio += rejections_met/random_walkers		
		end


		# Define the next α
		dE_dalpha = 2*(Eln-E*ln)
		α = α + 0.05#γ*dE_dalpha

		# Update the arrays for plotting
		energy[i] = E
		alpha[i] = α
		variance[i] = E2 - E^2
	end
	return energy,alpha,variance
end;

# ╔═╡ 20151696-5be9-4b68-afb7-8a7217c16402
energy,alpha,variance = iterate_alpha();

# ╔═╡ e5f66da9-ef79-4f82-9e1b-6ad7b27ec75d
energy

# ╔═╡ 64f872bd-cb8b-42de-8fa3-3270430e9a7f
alpha

# ╔═╡ 5ca1aa2b-099e-479c-88ad-7c9764fa6f85
begin
	plot(energy,
		label = false,
		title = "Energy expenctation Plot",
		linecolor = :blue,
		ylabel = "<E>",
		xlabel = "timestep"
	)
	
	hline!([π^2/8],
		label = false
	)
end

# ╔═╡ b21c39ff-68f7-4725-9bdc-4c180b5b368b
plot(variance,
	title = "Variance in energy",
	label = false,
	linecolor = :blue,
	ylabel = "Var(E)",
	xlabel = "timestep"
)

# ╔═╡ bb6cfb6d-325e-4913-9f77-3efcfcfc4368
plot(alpha,
	title = "Variational Parameter",
	label = false,
	linecolor = :blue,
	ylabel = "α",
	xlabel = "timestep"
)

# ╔═╡ 11dd75f3-e142-4712-914d-d8d594acada3
begin
	scatter(alpha,
		energy,
		title = "Energy vs α",
		xlabel = "α",
		ylabel = "<E>",
		label = false,
		markercolor = :blue
	)
	
	hline!([π^2/8],
		label = false
	)
end

# ╔═╡ 43a28afe-622c-4cfc-9679-9554994be6bf
variance

# ╔═╡ Cell order:
# ╠═7a0128c0-a13f-11eb-1098-ef9c00274a0d
# ╠═7c5e40c4-f84e-4527-a87e-c6ac395f5d9e
# ╠═007bb513-b173-4723-9d8b-2e858ae54dc6
# ╠═dcc6bdb1-daa8-44f1-a6ec-285ac1916179
# ╠═99eff10f-5b1e-42d9-9009-56e5aeb821c2
# ╠═cbfdf51f-721e-4172-b7d7-1911bf666163
# ╠═20151696-5be9-4b68-afb7-8a7217c16402
# ╠═e5f66da9-ef79-4f82-9e1b-6ad7b27ec75d
# ╠═64f872bd-cb8b-42de-8fa3-3270430e9a7f
# ╠═5ca1aa2b-099e-479c-88ad-7c9764fa6f85
# ╠═b21c39ff-68f7-4725-9bdc-4c180b5b368b
# ╠═bb6cfb6d-325e-4913-9f77-3efcfcfc4368
# ╠═11dd75f3-e142-4712-914d-d8d594acada3
# ╠═43a28afe-622c-4cfc-9679-9554994be6bf
