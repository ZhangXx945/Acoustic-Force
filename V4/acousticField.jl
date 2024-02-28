### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ 0802ac4a-b44b-11ee-1a5b-afdf3cf66b9b
"""
pressure field in (x,y,z) position in background medium
"""
function pressureField(x::Float64, y::Float64, z::Float64, ωb::Float64,
acoustic_model::FrequencySimulation, scattering_coefficient::Matrix{ComplexF64})
	ptc_No=length(acoustic_model.particles);
	dimension=typeof(acoustic_model.source.medium).parameters[2];
	cb=Float64(acoustic_model.source.medium.c);
	coef_order=Int(sqrt(length(scattering_coefficient[:,1])))-1;
	k=ωb/cb;
	pressure_field=0.0+0.0*im;
	for i in 1:ptc_No
		xx=x-acoustic_model.particles[i].shape.origin[1];
		yy=y-acoustic_model.particles[i].shape.origin[2];
		zz=z-acoustic_model.particles[i].shape.origin[3];
		r=sqrt(xx*xx+yy*yy+zz*zz);
		θ=atan(sqrt(xx*xx+yy*yy),zz);
		ϕ=atan(yy,xx);
		for n in 0:coef_order
			for m in -n:n
				pressure_field=pressure_field+scattering_coefficient[n*n+n+m+1,i]*hankel1(n,k*r)*ymn(n,m,θ,ϕ);
			end
		end
	end
	return pressure_field+acoustic_model.source.field([x,y,z],ωb);
end

# ╔═╡ 76254f0f-c39f-49b6-b5c3-7838b70776bd
"""
pressure field in (x,y,z) position inside particle ptc_which
"""
function pressureField(x::Float64, y::Float64, z::Float64, ωb::Float64,
acoustic_model::FrequencySimulation, scattering_coefficient::Matrix{ComplexF64},ptc_which::Int64)
	ptc_No=length(acoustic_model.particles);
	dimension=typeof(acoustic_model.source.medium).parameters[2];
	cp=Float64(acoustic_model.particles[ptc_which].medium.c);
	coef_order=Int(sqrt(length(scattering_coefficient[:,1])))-1;
	k=ωb/cp;
	pressure_field=0.0+0.0*im;
	for i in 1:ptc_No
		xx=x-acoustic_model.particles[i].shape.origin[1];
		yy=y-acoustic_model.particles[i].shape.origin[2];
		zz=z-acoustic_model.particles[i].shape.origin[3];
		r=sqrt(xx*xx+yy*yy+zz*zz);
		θ=atan(sqrt(xx*xx+yy*yy),zz);
		ϕ=atan(yy,xx);
		for n in 0:coef_order
			for m in -n:n
				if i==ptc_which
					pressure_field=pressure_field+scattering_coefficient[n*n+n+m+1,i]*hankel2(n,k*r)*ymn(n,m,θ,ϕ);
				else
					pressure_field=pressure_field+scattering_coefficient[n*n+n+m+1,i]*hankel1(n,k*r)*ymn(n,m,θ,ϕ);
				end
			end
		end
	end
	return pressure_field+acoustic_model.source.field([x,y,z],ωb);
end

# ╔═╡ 98847818-9ab1-4498-9931-adcad354c20f
"""
velocity field in (x,y,z) position in background medium
"""
function velocityField(x::Float64, y::Float64, z::Float64, ωb::Float64,
acoustic_model::FrequencySimulation, scattering_coefficient::Matrix{ComplexF64},δh::Float64)
	#=δh is suggessted to be 0.00005=#
	ρb=acoustic_model.source.medium.ρ;
	pxP1=pressureField(x+δh,y,z,ωb,acoustic_model,scattering_coefficient);#=the capital letter P in pxP1 denotes the positive step=#
	pxP2=pressureField(x+2*δh,y,z,ωb,acoustic_model,scattering_coefficient);
	pxN1=pressureField(x-δh,y,z,ωb,acoustic_model,scattering_coefficient);#=the capital letter P in pxP1 denotes the negative step=#
	pxN2=pressureField(x-2*δh,y,z,ωb,acoustic_model,scattering_coefficient);
	pyP1=pressureField(x,y+δh,z,ωb,acoustic_model,scattering_coefficient);
	pyP2=pressureField(x,y+2*δh,z,ωb,acoustic_model,scattering_coefficient);
	pyN1=pressureField(x,y-δh,z,ωb,acoustic_model,scattering_coefficient);
	pyN2=pressureField(x,y-2*δh,z,ωb,acoustic_model,scattering_coefficient);
	pzP1=pressureField(x,y,z+δh,ωb,acoustic_model,scattering_coefficient);
	pzP2=pressureField(x,y,z+2*δh,ωb,acoustic_model,scattering_coefficient);
	pzN1=pressureField(x,y,z-δh,ωb,acoustic_model,scattering_coefficient);
	pzN2=pressureField(x,y,z-2*δh,ωb,acoustic_model,scattering_coefficient);
	return -im/ρb/ωb*[(-pxP2+8*pxP1-8*pxN1+pxN2)/12/δh;(-pyP2+8*pyP1-8*pyN1+pyN2)/12/δh;(-pzP2+8*pzP1-8*pzN1+pzN2)/12/δh]
end

# ╔═╡ f8def69f-8085-42a1-8975-01facea50a7b
"""
velocity field in (x,y,z) position in background medium
"""
function velocityField(x::Float64, y::Float64, z::Float64, ωb::Float64,
acoustic_model::FrequencySimulation, scattering_coefficient::Matrix{ComplexF64},ptc_which::Int64,δh::Float64)
	#=δh is suggessted to be 0.00005=#
	ρp=acoustic_model.particles[ptc_which].medium.ρ;
	pxP1=pressureField(x+δh,y,z,ωb,acoustic_model,scattering_coefficient,ptc_which);#=the capital letter P in pxP1 denotes the positive step=#
	pxP2=pressureField(x+2*δh,y,z,ωb,acoustic_model,scattering_coefficient,ptc_which);
	pxN1=pressureField(x-δh,y,z,ωb,acoustic_model,scattering_coefficient,ptc_which);#=the capital letter P in pxP1 denotes the negative step=#
	pxN2=pressureField(x-2*δh,y,z,ωb,acoustic_model,scattering_coefficient,ptc_which);
	pyP1=pressureField(x,y+δh,z,ωb,acoustic_model,scattering_coefficient,ptc_which);
	pyP2=pressureField(x,y+2*δh,z,ωb,acoustic_model,scattering_coefficient,ptc_which);
	pyN1=pressureField(x,y-δh,z,ωb,acoustic_model,scattering_coefficient,ptc_which);
	pyN2=pressureField(x,y-2*δh,z,ωb,acoustic_model,scattering_coefficient,ptc_which);
	pzP1=pressureField(x,y,z+δh,ωb,acoustic_model,scattering_coefficient,ptc_which);
	pzP2=pressureField(x,y,z+2*δh,ωb,acoustic_model,scattering_coefficient,ptc_which);
	pzN1=pressureField(x,y,z-δh,ωb,acoustic_model,scattering_coefficient,ptc_which);
	pzN2=pressureField(x,y,z-2*δh,ωb,acoustic_model,scattering_coefficient,ptc_which);
	return -im/ρp/ωb*[(-pxP2+8*pxP1-8*pxN1+pxN2)/12/δh;(-pyP2+8*pyP1-8*pyN1+pyN2)/12/δh;(-pzP2+8*pzP1-8*pzN1+pzN2)/12/δh]
end

# ╔═╡ Cell order:
# ╠═0802ac4a-b44b-11ee-1a5b-afdf3cf66b9b
# ╠═76254f0f-c39f-49b6-b5c3-7838b70776bd
# ╠═98847818-9ab1-4498-9931-adcad354c20f
# ╠═f8def69f-8085-42a1-8975-01facea50a7b
