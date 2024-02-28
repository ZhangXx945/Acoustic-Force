### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ fee68378-cf1d-11ee-1b72-efa4b58f1637
"""
calculate force matrix
"""
function forceMatrix!(θsteps::Int64,ϕsteps::Int64,ωb::Float64,acoustic_model::FrequencySimulation,coef_order::Int64,δr::Float64,δh::Float64,δh2::Float64,force_matrix::Matrix{Float64})
	ptc_No=length(acoustic_model.particles);
	ptc_position=Matrix{Float64}(undef,ptc_No,3)
	for i in 1:ptc_No
		ptc_position[i,:].=acoustic_model.particles[i].shape.origin;
	end
	plus_δh2_position=deepcopy(ptc_position);
	plus_δh2_force=Matrix{Float64}(undef,ptc_No,3);
	plus_2δh2_position=deepcopy(ptc_position);
	plus_2δh2_force=Matrix{Float64}(undef,ptc_No,3);
	minus_δh2_position=deepcopy(ptc_position);
	minus_δh2_force=Matrix{Float64}(undef,ptc_No,3);
	minus_2δh2_position=deepcopy(ptc_position);
	minus_2δh2_force=Matrix{Float64}(undef,ptc_No,3);
	for i in 1:ptc_No
		for j in 1:2
			
			plus_δh2_position[i,j]=plus_δh2_position[i,j]+δh2;
			plus_δh2_model=buildModel2(acoustic_model,plus_δh2_position);
			plus_δh2_scattering_coefficient=getCoefficients(plus_δh2_model,ωb,coef_order);
			
			allParticleAcousticForce!(θsteps,ϕsteps,ωb,plus_δh2_model,plus_δh2_scattering_coefficient,δr,δh,plus_δh2_force);
			
			plus_2δh2_position[i,j]=plus_2δh2_position[i,j]+2*δh2;
			plus_2δh2_model=buildModel2(acoustic_model,plus_2δh2_position);
			plus_2δh2_scattering_coefficient=getCoefficients(plus_2δh2_model,ωb,coef_order);
			
			allParticleAcousticForce!(θsteps,ϕsteps,ωb,plus_2δh2_model,plus_2δh2_scattering_coefficient,δr,δh,plus_2δh2_force);
			
			minus_δh2_position[i,j]=minus_δh2_position[i,j]-δh2;
			minus_δh2_model=buildModel2(acoustic_model,minus_δh2_position);
			minus_δh2_scattering_coefficient=getCoefficients(minus_δh2_model,ωb,coef_order);
			
			allParticleAcousticForce!(θsteps,ϕsteps,ωb,minus_δh2_model,minus_δh2_scattering_coefficient,δr,δh,minus_δh2_force);
			
			minus_2δh2_position[i,j]=minus_2δh2_position[i,j]-2*δh2;
			minus_2δh2_model=buildModel2(acoustic_model,minus_2δh2_position);
			minus_2δh2_scattering_coefficient=getCoefficients(minus_2δh2_model,ωb,coef_order);
			
			allParticleAcousticForce!(θsteps,ϕsteps,ωb,minus_2δh2_model,minus_2δh2_scattering_coefficient,δr,δh,minus_2δh2_force);
			#====five points stencil to calculate force matrix element====#
			for m in 1:ptc_No
				for n in 1:2
					force_matrix[n+2*(m-1),j+2*(i-1)]=(-plus_2δh2_force[m,n]+8*plus_δh2_force[m,n]-8*minus_δh2_force[m,n]+minus_2δh2_force[m,n])/12/δh2;
				end
			end
			plus_δh2_position.=ptc_position;
			plus_2δh2_position.=ptc_position;
			minus_δh2_position.=ptc_position;
			minus_2δh2_position.=ptc_position;
		end
	end
end

# ╔═╡ Cell order:
# ╠═fee68378-cf1d-11ee-1b72-efa4b58f1637
