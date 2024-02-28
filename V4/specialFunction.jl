### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ 606cb93c-af5f-11ee-381b-e173f12ae223
#=spherical Hankel function of the first kind and the second kind, the spherical harmonics ymn=#
hankel1(n::Integer,x::Float64)=sf_bessel_jl(n,x)+im*sf_bessel_yl(n,x);

# ╔═╡ 2a9314b8-fb59-4dc1-8537-c5f0d2044d35
hankel2(n::Integer,x::Float64)=sf_bessel_jl(n,x)-im*sf_bessel_yl(n,x);

# ╔═╡ 4d3bd86d-eda0-40fa-a76f-fb208037673e
"""
Spherical harmonics
"""
function ymn(n::Integer,m::Integer,θ::Float64,ϕ::Float64) #spherical harmonics
if m>=0
return sf_legendre_sphPlm(n,m,cos(θ))*exp(im*m*ϕ);
elseif isodd(-m)
return -sf_legendre_sphPlm(n,-m,cos(θ))*exp(im*m*ϕ);
else
return sf_legendre_sphPlm(n,-m,cos(θ))*exp(im*m*ϕ);
end
end

# ╔═╡ Cell order:
# ╠═606cb93c-af5f-11ee-381b-e173f12ae223
# ╠═2a9314b8-fb59-4dc1-8537-c5f0d2044d35
# ╠═4d3bd86d-eda0-40fa-a76f-fb208037673e
