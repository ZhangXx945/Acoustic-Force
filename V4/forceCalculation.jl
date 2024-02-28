### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ 6c9cad10-b69c-11ee-0096-ad5fd458b0e4
"""
force density in (r+δr,θ,ϕ) point centered at the ptc_id particle
"""
function forceDensity(δr::Float64,θ::Float64,ϕ::Float64,ptc_id::Integer,ωb::Float64, acoustic_model::FrequencySimulation, scattering_coefficient::Matrix{ComplexF64},δh::Float64)
	R=acoustic_model.particles[ptc_id].shape.radius+δr;
	x0=acoustic_model.particles[ptc_id].shape.origin[1];
	y0=acoustic_model.particles[ptc_id].shape.origin[2];
	z0=acoustic_model.particles[ptc_id].shape.origin[3];
	sθ=sin(θ);cθ=cos(θ);sϕ=sin(ϕ);cϕ=cos(ϕ);
	x=x0+R*sθ*cϕ;
	y=y0+R*sθ*sϕ;
	z=z0+R*cθ;
	momentum_flux_density=Matrix{Float64}(undef,3,3);
	momentumFluxDensity!(x,y,z,ωb,acoustic_model,scattering_coefficient,δh,momentum_flux_density);
	force_density1=(momentum_flux_density[1,1]*sθ*cϕ+momentum_flux_density[1,2]*sθ*sϕ+momentum_flux_density[1,3]*cθ)*R*R*sθ;
	force_density2=(momentum_flux_density[2,1]*sθ*cϕ+momentum_flux_density[2,2]*sθ*sϕ+momentum_flux_density[2,3]*cθ)*R*R*sθ;
	force_density3=(momentum_flux_density[3,1]*sθ*cϕ+momentum_flux_density[3,2]*sθ*sϕ+momentum_flux_density[3,3]*cθ)*R*R*sθ;
	momentum_flux_density=Nothing;
	return [force_density1 force_density2 force_density3]
end

# ╔═╡ aae0324b-391b-4a49-a073-f82c7d7235b7
"""
parallel computation on part of closed surface integral of force density
"""
function forcePart1(ptc_id::Int64,θsteps::Int64,ϕsteps::Int64,θ_lower_limit::Float64,θ_upper_limit::Float64,ϕ_lower_limit::Float64,ϕ_upper_limit::Float64,ωb::Float64,acoustic_model::FrequencySimulation,scattering_coefficient::Matrix{ComplexF64},δr::Float64,δh::Float64)
	xforce=0.0;yforce=0.0;zforce=0.0;
	δθ=(θ_upper_limit-θ_lower_limit)/θsteps;
	δϕ=(ϕ_upper_limit-ϕ_lower_limit)/ϕsteps;
	for i in 1:θsteps
		θ=θ_lower_limit+(i-1)*δθ;
		for j in 1:ϕsteps
			ϕ=ϕ_lower_limit+(j-1)*δϕ;
			force_density=forceDensity(δr,θ,ϕ,ptc_id,ωb,acoustic_model,scattering_coefficient,δh);
			xforce=xforce+force_density[1];
			yforce=yforce+force_density[2];
			zforce=zforce+force_density[3];
		end
	end
	return -δθ*δϕ*[xforce yforce zforce]
end

# ╔═╡ 83c0d01b-3969-422b-a6b3-43434b247972
"""
parallel computation on part of closed surface integral of force density
"""
function forcePart2(ptc_id::Int64,θsteps::Int64,ϕsteps::Int64,θ_lower_limit::Float64,θ_upper_limit::Float64,ϕ_lower_limit::Float64,ϕ_upper_limit::Float64,ωb::Float64,acoustic_model::FrequencySimulation,scattering_coefficient::Matrix{ComplexF64},δr::Float64,δh::Float64)
	xforce=0.0;yforce=0.0;zforce=0.0;
	δθ=(θ_upper_limit-θ_lower_limit)/θsteps;
	δϕ=(ϕ_upper_limit-ϕ_lower_limit)/ϕsteps;
	for i in 1:θsteps
		θ=θ_lower_limit+(i-1)*δθ;
		for j in 1:ϕsteps
			ϕ=ϕ_lower_limit+(j-1)*δϕ;
			force_density=forceDensity(δr,θ,ϕ,ptc_id,ωb,acoustic_model,scattering_coefficient,δh);
			xforce=xforce+force_density[1];
			yforce=yforce+force_density[2];
			zforce=zforce+force_density[3];
		end
	end
	return -δθ*δϕ*[xforce yforce zforce]
end

# ╔═╡ 82d82a2b-92a3-4610-990e-c9d1933cb6d0
"""
parallel computation on part of closed surface integral of force density
"""
function forcePart3(ptc_id::Int64,θsteps::Int64,ϕsteps::Int64,θ_lower_limit::Float64,θ_upper_limit::Float64,ϕ_lower_limit::Float64,ϕ_upper_limit::Float64,ωb::Float64,acoustic_model::FrequencySimulation,scattering_coefficient::Matrix{ComplexF64},δr::Float64,δh::Float64)
	xforce=0.0;yforce=0.0;zforce=0.0;
	δθ=(θ_upper_limit-θ_lower_limit)/θsteps;
	δϕ=(ϕ_upper_limit-ϕ_lower_limit)/ϕsteps;
	for i in 1:θsteps
		θ=θ_lower_limit+(i-1)*δθ;
		for j in 1:ϕsteps
			ϕ=ϕ_lower_limit+(j-1)*δϕ;
			force_density=forceDensity(δr,θ,ϕ,ptc_id,ωb,acoustic_model,scattering_coefficient,δh);
			xforce=xforce+force_density[1];
			yforce=yforce+force_density[2];
			zforce=zforce+force_density[3];
		end
	end
	return -δθ*δϕ*[xforce yforce zforce]
end

# ╔═╡ 36e6de66-041c-4e8a-8d99-0d822d1e212a
"""
parallel computation on part of closed surface integral of force density
"""
function forcePart4(ptc_id::Int64,θsteps::Int64,ϕsteps::Int64,θ_lower_limit::Float64,θ_upper_limit::Float64,ϕ_lower_limit::Float64,ϕ_upper_limit::Float64,ωb::Float64,acoustic_model::FrequencySimulation,scattering_coefficient::Matrix{ComplexF64},δr::Float64,δh::Float64)
	xforce=0.0;yforce=0.0;zforce=0.0;
	δθ=(θ_upper_limit-θ_lower_limit)/θsteps;
	δϕ=(ϕ_upper_limit-ϕ_lower_limit)/ϕsteps;
	for i in 1:θsteps
		θ=θ_lower_limit+(i-1)*δθ;
		for j in 1:ϕsteps
			ϕ=ϕ_lower_limit+(j-1)*δϕ;
			force_density=forceDensity(δr,θ,ϕ,ptc_id,ωb,acoustic_model,scattering_coefficient,δh);
			xforce=xforce+force_density[1];
			yforce=yforce+force_density[2];
			zforce=zforce+force_density[3];
		end
	end
	return -δθ*δϕ*[xforce yforce zforce]
end

# ╔═╡ 5de24b35-e540-43db-9571-58572407471b
"""
parallel computation on part of closed surface integral of force density
"""
function forcePart5(ptc_id::Int64,θsteps::Int64,ϕsteps::Int64,θ_lower_limit::Float64,θ_upper_limit::Float64,ϕ_lower_limit::Float64,ϕ_upper_limit::Float64,ωb::Float64,acoustic_model::FrequencySimulation,scattering_coefficient::Matrix{ComplexF64},δr::Float64,δh::Float64)
	xforce=0.0;yforce=0.0;zforce=0.0;
	δθ=(θ_upper_limit-θ_lower_limit)/θsteps;
	δϕ=(ϕ_upper_limit-ϕ_lower_limit)/ϕsteps;
	for i in 1:θsteps
		θ=θ_lower_limit+(i-1)*δθ;
		for j in 1:ϕsteps
			ϕ=ϕ_lower_limit+(j-1)*δϕ;
			force_density=forceDensity(δr,θ,ϕ,ptc_id,ωb,acoustic_model,scattering_coefficient,δh);
			xforce=xforce+force_density[1];
			yforce=yforce+force_density[2];
			zforce=zforce+force_density[3];
		end
	end
	return -δθ*δϕ*[xforce yforce zforce]
end

# ╔═╡ 53122b44-7fc7-46a1-b9fa-77c5036f38bc
"""
parallel computation on part of closed surface integral of force density
"""
function forcePart6(ptc_id::Int64,θsteps::Int64,ϕsteps::Int64,θ_lower_limit::Float64,θ_upper_limit::Float64,ϕ_lower_limit::Float64,ϕ_upper_limit::Float64,ωb::Float64,acoustic_model::FrequencySimulation,scattering_coefficient::Matrix{ComplexF64},δr::Float64,δh::Float64)
	xforce=0.0;yforce=0.0;zforce=0.0;
	δθ=(θ_upper_limit-θ_lower_limit)/θsteps;
	δϕ=(ϕ_upper_limit-ϕ_lower_limit)/ϕsteps;
	for i in 1:θsteps
		θ=θ_lower_limit+(i-1)*δθ;
		for j in 1:ϕsteps
			ϕ=ϕ_lower_limit+(j-1)*δϕ;
			force_density=forceDensity(δr,θ,ϕ,ptc_id,ωb,acoustic_model,scattering_coefficient,δh);
			xforce=xforce+force_density[1];
			yforce=yforce+force_density[2];
			zforce=zforce+force_density[3];
		end
	end
	return -δθ*δϕ*[xforce yforce zforce]
end

# ╔═╡ d143e569-46d6-481c-9d73-b6040e8222a7
"""
parallel computation on part of closed surface integral of force density
"""
function forcePart7(ptc_id::Int64,θsteps::Int64,ϕsteps::Int64,θ_lower_limit::Float64,θ_upper_limit::Float64,ϕ_lower_limit::Float64,ϕ_upper_limit::Float64,ωb::Float64,acoustic_model::FrequencySimulation,scattering_coefficient::Matrix{ComplexF64},δr::Float64,δh::Float64)
	xforce=0.0;yforce=0.0;zforce=0.0;
	δθ=(θ_upper_limit-θ_lower_limit)/θsteps;
	δϕ=(ϕ_upper_limit-ϕ_lower_limit)/ϕsteps;
	for i in 1:θsteps
		θ=θ_lower_limit+(i-1)*δθ;
		for j in 1:ϕsteps
			ϕ=ϕ_lower_limit+(j-1)*δϕ;
			force_density=forceDensity(δr,θ,ϕ,ptc_id,ωb,acoustic_model,scattering_coefficient,δh);
			xforce=xforce+force_density[1];
			yforce=yforce+force_density[2];
			zforce=zforce+force_density[3];
		end
	end
	return -δθ*δϕ*[xforce yforce zforce]
end

# ╔═╡ 74fcf743-b96e-4764-bd6a-979684e59b88
"""
parallel computation on part of closed surface integral of force density
"""
function forcePart8(ptc_id::Int64,θsteps::Int64,ϕsteps::Int64,θ_lower_limit::Float64,θ_upper_limit::Float64,ϕ_lower_limit::Float64,ϕ_upper_limit::Float64,ωb::Float64,acoustic_model::FrequencySimulation,scattering_coefficient::Matrix{ComplexF64},δr::Float64,δh::Float64)
	xforce=0.0;yforce=0.0;zforce=0.0;
	δθ=(θ_upper_limit-θ_lower_limit)/θsteps;
	δϕ=(ϕ_upper_limit-ϕ_lower_limit)/ϕsteps;
	for i in 1:θsteps
		θ=θ_lower_limit+(i-1)*δθ;
		for j in 1:ϕsteps
			ϕ=ϕ_lower_limit+(j-1)*δϕ;
			force_density=forceDensity(δr,θ,ϕ,ptc_id,ωb,acoustic_model,scattering_coefficient,δh);
			xforce=xforce+force_density[1];
			yforce=yforce+force_density[2];
			zforce=zforce+force_density[3];
		end
	end
	return -δθ*δϕ*[xforce yforce zforce]
end

# ╔═╡ 1a93f19d-8546-49a1-ae60-bc9a4170a7ca
"""
single-threaded routine to calculate acoutic force
"""
function forceT1(ptc_id::Int64,θsteps::Int64,ϕsteps::Int64,ωb::Float64,acoustic_model::FrequencySimulation,scattering_coefficient::Matrix{ComplexF64},δr::Float64,δh::Float64)
return forcePart1(ptc_id,θsteps,ϕsteps,0.0,1π,0.0,2π,ωb,acoustic_model,scattering_coefficient,δr,δh)
end

# ╔═╡ 2aaae57c-b9dc-4a6b-9ac2-98a60d03d18e
"""
2-threaded routine to calculate acoutic force
"""
function forceT2(ptc_id::Int64,θsteps::Int64,ϕsteps::Int64,ωb::Float64,acoustic_model::FrequencySimulation,scattering_coefficient::Matrix{ComplexF64},δr::Float64,δh::Float64)
	ϕsteps_tmp=Int(trunc(ϕsteps/2));
	force_part_1=Threads.@spawn forcePart1(ptc_id,θsteps,ϕsteps_tmp,0.0,1π,0.0,1π,ωb,acoustic_model,scattering_coefficient,δr,δh);
	force_part_2=Threads.@spawn forcePart2(ptc_id,θsteps,ϕsteps_tmp,0.0,1π,1π,2π,ωb,acoustic_model,scattering_coefficient,δr,δh);
	return fetch(force_part_1) + fetch(force_part_2)
end

# ╔═╡ 3f38934f-a010-4de4-a6d4-87a0f3125790
"""
4-threaded routine to calculate acoutic force
"""
function forceT4(ptc_id::Int64,θsteps::Int64,ϕsteps::Int64,ωb::Float64,acoustic_model::FrequencySimulation,scattering_coefficient::Matrix{ComplexF64},δr::Float64,δh::Float64)
	θsteps_tmp=Int(trunc(θsteps/2));
	ϕsteps_tmp=Int(trunc(ϕsteps/2));
	force_part_1=Threads.@spawn forcePart1(ptc_id,θsteps_tmp,ϕsteps_tmp,0.0,0.5π,0.0,1π,ωb,acoustic_model,scattering_coefficient,δr,δh);
	force_part_2=Threads.@spawn forcePart2(ptc_id,θsteps_tmp,ϕsteps_tmp,0.0,0.5π,1π,2π,ωb,acoustic_model,scattering_coefficient,δr,δh);
	force_part_3=Threads.@spawn forcePart3(ptc_id,θsteps_tmp,ϕsteps_tmp,0.5π,1π,0.0,1π,ωb,acoustic_model,scattering_coefficient,δr,δh);
	force_part_4=Threads.@spawn forcePart4(ptc_id,θsteps_tmp,ϕsteps_tmp,0.5π,1π,1π,2π,ωb,acoustic_model,scattering_coefficient,δr,δh);
	return fetch(force_part_1) + fetch(force_part_2)+fetch(force_part_3)+fetch(force_part_4)
end

# ╔═╡ 736cfa08-8767-409c-8653-7c5ad0cb35d0
"""
6-threaded routine to calculate acoutic force
"""
function forceT6(ptc_id::Int64,θsteps::Int64,ϕsteps::Int64,ωb::Float64,acoustic_model::FrequencySimulation,scattering_coefficient::Matrix{ComplexF64},δr::Float64,δh::Float64)
	θsteps_tmp=Int(trunc(θsteps/2));
	ϕsteps_tmp=Int(trunc(ϕsteps/3));
	force_part_1=Threads.@spawn forcePart1(ptc_id,θsteps_tmp,ϕsteps_tmp,0.0,0.5π,0.0,2π/3,ωb,acoustic_model,scattering_coefficient,δr,δh);
	force_part_2=Threads.@spawn forcePart2(ptc_id,θsteps_tmp,ϕsteps_tmp,0.0,0.5π,2π/3,4π/3,ωb,acoustic_model,scattering_coefficient,δr,δh);
	force_part_3=Threads.@spawn forcePart3(ptc_id,θsteps_tmp,ϕsteps_tmp,0.0,0.5π,4π/3,2π,ωb,acoustic_model,scattering_coefficient,δr,δh);
	force_part_4=Threads.@spawn forcePart4(ptc_id,θsteps_tmp,ϕsteps_tmp,0.5π,1π,0.0,2π/3,ωb,acoustic_model,scattering_coefficient,δr,δh);
	force_part_5=Threads.@spawn forcePart5(ptc_id,θsteps_tmp,ϕsteps_tmp,0.5π,1π,2π/3,4π/3,ωb,acoustic_model,scattering_coefficient,δr,δh);
	force_part_6=Threads.@spawn forcePart6(ptc_id,θsteps_tmp,ϕsteps_tmp,0.5π,1π,4π/3,2π,ωb,acoustic_model,scattering_coefficient,δr,δh);
	return fetch(force_part_1) + fetch(force_part_2)+fetch(force_part_3)+fetch(force_part_4)+fetch(force_part_5)+fetch(force_part_6)
end

# ╔═╡ 563c23e6-6897-44da-acc4-8e0c9ed3bc02
"""
8-threaded routine to calculate acoutic force
"""
function forceT8(ptc_id::Int64,θsteps::Int64,ϕsteps::Int64,ωb::Float64,acoustic_model::FrequencySimulation,scattering_coefficient::Matrix{ComplexF64},δr::Float64,δh::Float64)
	θsteps_tmp=Int(trunc(θsteps/2));
	ϕsteps_tmp=Int(trunc(ϕsteps/4));
	force_part_1=Threads.@spawn forcePart1(ptc_id,θsteps_tmp,ϕsteps_tmp,0.0,0.5π,0.0,π/2,ωb,acoustic_model,scattering_coefficient,δr,δh);
	force_part_2=Threads.@spawn forcePart2(ptc_id,θsteps_tmp,ϕsteps_tmp,0.0,0.5π,π/2,1π,ωb,acoustic_model,scattering_coefficient,δr,δh);
	force_part_3=Threads.@spawn forcePart3(ptc_id,θsteps_tmp,ϕsteps_tmp,0.0,0.5π,1π,3π/2,ωb,acoustic_model,scattering_coefficient,δr,δh);
	force_part_4=Threads.@spawn forcePart4(ptc_id,θsteps_tmp,ϕsteps_tmp,0.0,0.5π,3π/2,2π,ωb,acoustic_model,scattering_coefficient,δr,δh);
	force_part_5=Threads.@spawn forcePart5(ptc_id,θsteps_tmp,ϕsteps_tmp,0.5π,1π,0.0,π/2,ωb,acoustic_model,scattering_coefficient,δr,δh);
	force_part_6=Threads.@spawn forcePart6(ptc_id,θsteps_tmp,ϕsteps_tmp,0.5π,1π,π/2,1π,ωb,acoustic_model,scattering_coefficient,δr,δh);
	force_part_7=Threads.@spawn forcePart7(ptc_id,θsteps_tmp,ϕsteps_tmp,0.5π,1π,1π,3π/2,ωb,acoustic_model,scattering_coefficient,δr,δh);
	force_part_8=Threads.@spawn forcePart8(ptc_id,θsteps_tmp,ϕsteps_tmp,0.5π,1π,3π/2,2π,ωb,acoustic_model,scattering_coefficient,δr,δh);
	return fetch(force_part_1) + fetch(force_part_2)+fetch(force_part_3)+fetch(force_part_4)+fetch(force_part_5)+fetch(force_part_6)+fetch(force_part_7)+fetch(force_part_8)
end

# ╔═╡ 2f291074-d932-40ce-8b74-473e03626274
"""
acoustic force on a single particle ptc_id
"""
function acousticForce(ptc_id::Int64,θsteps::Int64,ϕsteps::Int64,ωb::Float64,acoustic_model::FrequencySimulation,scattering_coefficient::Matrix{ComplexF64},δr::Float64,δh::Float64)
	num_threads = Threads.nthreads();
	if num_threads==1
		println("Single-threaded computation")
		return forceT1(ptc_id,θsteps,ϕsteps,ωb,acoustic_model,scattering_coefficient,δr,δh)
	elseif 2<=num_threads&&num_threads<4
		println("2-threaded parallel computation")
		return forceT2(ptc_id,θsteps,ϕsteps,ωb,acoustic_model,scattering_coefficient,δr,δh)
	elseif 4<=num_threads&&num_threads<6
		println("4-threaded parallel computation")
		return forceT4(ptc_id,θsteps,ϕsteps,ωb,acoustic_model,scattering_coefficient,δr,δh)
	elseif 6<=num_threads&&num_threads<8
		println("6-threaded parallel computation")
		return forceT6(ptc_id,θsteps,ϕsteps,ωb,acoustic_model,scattering_coefficient,δr,δh)
	else
	
		println("8-threaded parallel computation")
		return forceT8(ptc_id,θsteps,ϕsteps,ωb,acoustic_model,scattering_coefficient,δr,δh)
	end
end

# ╔═╡ e0d47cc6-565c-4f97-a2ad-9004a6b6dcfd
function allParticleAcousticForce!(θsteps::Int64,ϕsteps::Int64,ωb::Float64,acoustic_model::FrequencySimulation,scattering_coefficient::Matrix{ComplexF64},δr::Float64,δh::Float64,all_particle_force::Matrix{Float64})
	ptc_No=length(acoustic_model.particles);
	for i in 1:ptc_No
		all_particle_force[i,:]=acousticForce(i,θsteps,ϕsteps,ωb,acoustic_model,scattering_coefficient,δr,δh);
	end
end

# ╔═╡ Cell order:
# ╠═6c9cad10-b69c-11ee-0096-ad5fd458b0e4
# ╠═aae0324b-391b-4a49-a073-f82c7d7235b7
# ╠═83c0d01b-3969-422b-a6b3-43434b247972
# ╠═82d82a2b-92a3-4610-990e-c9d1933cb6d0
# ╠═36e6de66-041c-4e8a-8d99-0d822d1e212a
# ╠═5de24b35-e540-43db-9571-58572407471b
# ╠═53122b44-7fc7-46a1-b9fa-77c5036f38bc
# ╠═d143e569-46d6-481c-9d73-b6040e8222a7
# ╠═74fcf743-b96e-4764-bd6a-979684e59b88
# ╠═1a93f19d-8546-49a1-ae60-bc9a4170a7ca
# ╠═2aaae57c-b9dc-4a6b-9ac2-98a60d03d18e
# ╠═3f38934f-a010-4de4-a6d4-87a0f3125790
# ╠═736cfa08-8767-409c-8653-7c5ad0cb35d0
# ╠═563c23e6-6897-44da-acc4-8e0c9ed3bc02
# ╠═2f291074-d932-40ce-8b74-473e03626274
# ╠═e0d47cc6-565c-4f97-a2ad-9004a6b6dcfd
