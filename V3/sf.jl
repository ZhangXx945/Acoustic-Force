### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ ed3254e0-7b23-11ee-3e45-7dfc3a3015bf
hk(n::Integer,x::Float64)=sf_bessel_jl(n,x)+im*sf_bessel_yl(n,x);

# ╔═╡ 653b4957-8342-4039-80f0-9276ada454b1
function ymn(n::Integer,m::Integer,θ::Float64,ϕ::Float64) 	#spherical harmonics
if m>=0
	return sf_legendre_sphPlm(n,m,cos(θ))*exp(im*m*ϕ);
elseif isodd(-m)
	return -sf_legendre_sphPlm(n,-m,cos(θ))*exp(im*m*ϕ);
else
	return sf_legendre_sphPlm(n,-m,cos(θ))*exp(im*m*ϕ);
end
end

# ╔═╡ Cell order:
# ╠═ed3254e0-7b23-11ee-3e45-7dfc3a3015bf
# ╠═653b4957-8342-4039-80f0-9276ada454b1