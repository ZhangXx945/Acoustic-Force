### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ c1a76103-2c07-4f4e-824b-9f5367bfe833
"""
build the acoustic model
"""
function buildModel(wave_No::Int64,ρb::Float64,cb::Float64,inc_amplitude::Array{Float64},inc_direction::Matrix{Float64},inc_position::Matrix{Float64},ptc_No::Int64,radius::Array{Float64},ρp::Array{Float64},cp::Array{Float64},ptc_position::Matrix{Float64})
	#=set incident waves and background medium=#
	bg_property=Acoustic(3; ρ = ρb, c = cb);
	waves=plane_source(bg_property;amplitude=inc_amplitude[1],direction=inc_direction[1,:],position = inc_position[1,:]);
	for i in 2:wave_No
		waves=waves+plane_source(bg_property;amplitude=inc_amplitude[i],direction=inc_direction[i,:],position = inc_position[i,:]);
	end
	#=set particle property=#
	particles=[Particle(Acoustic(3; ρ = ρp[1], c = cp[1]),Sphere(ptc_position[1,:],radius[1]))];
	for i in 2:ptc_No
		push!(particles,Particle(Acoustic(3; ρ = ρp[i], c = cp[i]),Sphere(ptc_position[i,:],radius[i])));
	end
	acoustic_model=FrequencySimulation(particles,waves);
	particles=Nothing;
	return acoustic_model
end

# ╔═╡ f965d049-4e3e-406c-becb-4959bf695629
function buildModel2(acoustic_model::FrequencySimulation,ptc_position::Matrix{Float64})
	ptc_No=length(acoustic_model.particles);
	particles=[Particle(Acoustic(3;ρ=acoustic_model.particles[1].medium.ρ,c=acoustic_model.particles[1].medium.c),Sphere(ptc_position[1,:],acoustic_model.particles[1].shape.radius))];
	for i in 2:ptc_No
		push!(particles,Particle(Acoustic(3;ρ=acoustic_model.particles[i].medium.ρ,c=acoustic_model.particles[i].medium.c),Sphere(ptc_position[i,:],acoustic_model.particles[i].shape.radius)))
	end
	return FrequencySimulation(particles,acoustic_model.source)
end

# ╔═╡ 590dc1b6-ae2d-4eb2-b19f-c81db5e9d3f5
"""
get scattering coefficients
"""
function getCoefficients(acoustic_model::FrequencySimulation,ωb::Float64, coef_order::Int64)
	return basis_coefficients(acoustic_model,ωb,basis_order=coef_order)
end

# ╔═╡ Cell order:
# ╠═c1a76103-2c07-4f4e-824b-9f5367bfe833
# ╠═f965d049-4e3e-406c-becb-4959bf695629
# ╠═590dc1b6-ae2d-4eb2-b19f-c81db5e9d3f5
