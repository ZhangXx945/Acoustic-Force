### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ d6e67220-b5e4-11ee-29d5-4f08d9c7653e
"""
acoustic momentum flux density in position (x,y,z) in background medium
stress_tensor is:
Txx Txy Txz
Tyx Tyy Tyz
Tzx Tzy Tzz
"""
function momentumFluxDensity!(x::Float64, y::Float64, z::Float64, ωb::Float64,
acoustic_model::FrequencySimulation, scattering_coefficient::Matrix{ComplexF64},δh::Float64,momentum_flux_density::Matrix{Float64})
	ρb=acoustic_model.source.medium.ρ;
	cb=acoustic_model.source.medium.c;
	v=velocityField(x,y,z,ωb,acoustic_model,scattering_coefficient,δh);
	cjv=conj.(v);
	p=pressureField(x,y,z,ωb,acoustic_model,scattering_coefficient);
	cjp=conj(p);
	v_2=real(sum(v.*cjv));
	p_2_coef=real(p*cjp/2/ρb/cb/cb);
	momentum_flux_density[1,1]=0.5*(ρb*(real(v[1]*cjv[1])-0.5*v_2)+p_2_coef);
	momentum_flux_density[1,2]=0.5*ρb*real(v[1]*cjv[2]);
	momentum_flux_density[1,3]=0.5*ρb*real(v[1]*cjv[3]);
	momentum_flux_density[2,1]=0.5*ρb*real(v[2]*cjv[1]);
	momentum_flux_density[2,2]=0.5*(ρb*(real(v[2]*cjv[2])-0.5*v_2)+p_2_coef);
	momentum_flux_density[2,3]=0.5*ρb*real(v[2]*cjv[3]);
	momentum_flux_density[3,1]=0.5*ρb*real(v[3]*cjv[1]);
	momentum_flux_density[3,2]=0.5*ρb*real(v[3]*cjv[2]);
	momentum_flux_density[3,3]=0.5*(ρb*(real(v[3]*cjv[3])-0.5*v_2)+p_2_coef);
end

# ╔═╡ 1d31b6ed-8348-435d-b6b1-b7be05fc2afd
"""
acoustic momentum flux density in position (x,y,z) inside particle ptc_which
stress_tensor is:
Txx Txy Txz
Tyx Tyy Tyz
Tzx Tzy Tzz
"""
function momentumFluxDensity!(x::Float64, y::Float64, z::Float64, ωb::Float64,
acoustic_model::FrequencySimulation, scattering_coefficient::Matrix{ComplexF64},ptc_which::Int64,δh::Float64,momentum_flux_density::Matrix{Float64})
	ρp=acoustic_model.particles[ptc_which].medium.ρ;
	cp=acoustic_model.particles[ptc_which].medium.c;
	v=velocityField(x,y,z,ωb,acoustic_model,scattering_coefficient,ptc_which,δh);
	cjv=conj.(v);
	p=pressureField(x,y,z,ωb,acoustic_model,scattering_coefficient,ptc_which);
	cjp=conj(p);
	v_2=real(sum(v.*cjv));
	p_2_coef=real(p*cjp/2/ρp/cp/cp);
	momentum_flux_density[1,1]=0.5*(ρp*(real(v[1]*cjv[1])-0.5*v_2)+p_2_coef);
	momentum_flux_density[1,2]=0.5*ρp*real(v[1]*cjv[2]);
	momentum_flux_density[1,3]=0.5*ρp*real(v[1]*cjv[3]);
	momentum_flux_density[2,1]=0.5*ρp*real(v[2]*cjv[1]);
	momentum_flux_density[2,2]=0.5*(ρp*(real(v[2]*cjv[2])-0.5*v_2)+p_2_coef);
	momentum_flux_density[2,3]=0.5*ρp*real(v[2]*cjv[3]);
	momentum_flux_density[3,1]=0.5*ρp*real(v[3]*cjv[1]);
	momentum_flux_density[3,2]=0.5*ρp*real(v[3]*cjv[2]);
	momentum_flux_density[3,3]=0.5*(ρp*(real(v[3]*cjv[3])-0.5*v_2)+p_2_coef);
end

# ╔═╡ ee7d9514-5cb6-4fd7-a950-38ca7f62777e
"""
acoustic energy density in position (x,y,z)
"""
function energyDensity(x::Float64, y::Float64, z::Float64, ωb::Float64,
acoustic_model::FrequencySimulation, scattering_coefficient::Matrix{ComplexF64},δh::Float64)
	ρb=acoustic_model.source.medium.ρ;
	cb=acoustic_model.source.medium.c;
	v=velocityField(x,y,z,ωb,acoustic_model,scattering_coefficient,δh);
	cjv=conj.(v);
	p=pressureField(x,y,z,ωb,acoustic_model,scattering_coefficient);
	cjp=conj(p);
	β=1/ρb/cb/cb;
	return 0.25*(β*p*cjp+ρb*sum(v.*cjv))
end

# ╔═╡ 7dcabb5a-f488-4eab-808d-41aff0e45928
function energyDensity(x::Float64, y::Float64, z::Float64, ωb::Float64,
acoustic_model::FrequencySimulation, scattering_coefficient::Matrix{ComplexF64},ptc_which::Int64,δh::Float64)
	ρp=acoustic_model.particles[ptc_which].medium.ρ;
	cp=acoustic_model.particles[ptc_which].medium.c;
	v=velocityField(x,y,z,ωb,acoustic_model,scattering_coefficient,ptc_which,δh);
	cjv=conj.(v);
	p=pressureField(x,y,z,ωb,acoustic_model,scattering_coefficient,ptc_which);
	cjp=conj(p);
	β=1/ρp/cp/cp;
	return 0.25*(β*p*cjp+ρb*sum(v.*cjv))
end

# ╔═╡ e62ff277-8ce7-4694-8467-0d2a3fa9e89b
"""
acoustic energy flux density in position (x,y,z)
"""
function energyFluxDensity!(x::Float64, y::Float64, z::Float64, ωb::Float64,
acoustic_model::FrequencySimulation, scattering_coefficient::Matrix{ComplexF64},δh::Float64,energy_flux_density::Array{Float64})
	cjp=conj(pressureField(x,y,z,ωb,acoustic_model,scattering_coefficient));
	v=velocityField(x,y,z,ωb,acoustic_model,scattering_coefficient,δh);
	energy_flux_density[1]=0.5*real(cjp*v[1]);
	energy_flux_density[2]=0.5*real(cjp*v[2]);
	energy_flux_density[3]=0.5*real(cjp*v[3]);
end

# ╔═╡ d42f8b3b-4ea3-4837-b046-8ec2d52c705d
function energyFluxDensity!(x::Float64, y::Float64, z::Float64, ωb::Float64,
acoustic_model::FrequencySimulation, scattering_coefficient::Matrix{ComplexF64},ptc_which::Int64,δh::Float64,energy_flux_density::Array{Float64})
	cjp=conj(pressureField(x,y,z,ωb,acoustic_model,scattering_coefficient,ptc_which));
	v=velocityField(x,y,z,ωb,acoustic_model,scattering_coefficient,ptc_which,δh);
	energy_flux_density[1]=0.5*real(cjp*v[1]);
	energy_flux_density[2]=0.5*real(cjp*v[2]);
	energy_flux_density[3]=0.5*real(cjp*v[3]);
end

# ╔═╡ 383a3579-88bb-44d1-a32f-3e8b5adcc0a8
"""
acoustic spin flux density in position (x,y,z)
"""
function spinDensity!(x::Float64, y::Float64, z::Float64, ωb::Float64,
acoustic_model::FrequencySimulation, scattering_coefficient::Matrix{ComplexF64},δh::Float64,spin_density::Array{Float64})
	ρb=acoustic_model.source.medium.ρ;
	v=velocityField(x,y,z,ωb,acoustic_model,scattering_coefficient,δh);
	cjv=conj.(v);
	spin_density[1]=ρb/2/ωb*imag(cjv[2]*v[3]-cjv[3]*v[2]);
	spin_density[2]=ρb/2/ωb*imag(cjv[3]*v[1]-cjv[1]*v[3]);
	spin_density[3]=ρb/2/ωb*imag(cjv[1]*v[2]-cjv[2]*v[1]);
end

# ╔═╡ f09ed878-d870-4cae-acde-e62b30f1960f
function spinDensity!(x::Float64, y::Float64, z::Float64, ωb::Float64,
acoustic_model::FrequencySimulation, scattering_coefficient::Matrix{ComplexF64},ptc_which::Int64,δh::Float64,spin_density::Array{Float64})
	ρp=acoustic_model.particles[ptc_which].medium.ρ;
	v=velocityField(x,y,z,ωb,acoustic_model,scattering_coefficient,ptc_which,δh);
	cjv=conj.(v);
	spin_density[1]=ρp/2/ωb*imag(cjv[2]*v[3]-cjv[3]*v[2]);
	spin_density[2]=ρp/2/ωb*imag(cjv[3]*v[1]-cjv[1]*v[3]);
	spin_density[3]=ρp/2/ωb*imag(cjv[1]*v[2]-cjv[2]*v[1]);
end

# ╔═╡ Cell order:
# ╠═d6e67220-b5e4-11ee-29d5-4f08d9c7653e
# ╠═1d31b6ed-8348-435d-b6b1-b7be05fc2afd
# ╠═ee7d9514-5cb6-4fd7-a950-38ca7f62777e
# ╠═7dcabb5a-f488-4eab-808d-41aff0e45928
# ╠═e62ff277-8ce7-4694-8467-0d2a3fa9e89b
# ╠═d42f8b3b-4ea3-4837-b046-8ec2d52c705d
# ╠═383a3579-88bb-44d1-a32f-3e8b5adcc0a8
# ╠═f09ed878-d870-4cae-acde-e62b30f1960f
