### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 41b86c42-65a7-11ee-0959-ddd4537ec9be
using GSL

# ╔═╡ 24689bb2-8862-442f-9b82-03b2270c3ed6
using MultipleScattering

# ╔═╡ 4e14a103-6c4f-4bb0-8acc-c0be754be0c7
using Plots

# ╔═╡ 57493ea3-4975-4ddf-8955-06ecbb1a03ff
using PyPlot

# ╔═╡ 5bb8c76a-556b-4f2d-9af5-4c806220cba5
using LinearAlgebra

# ╔═╡ ea5429e2-5e91-4003-bc4e-5ee1f1980056
using DelimitedFiles

# ╔═╡ e245f8bd-f235-4ec6-b4e2-e405d650bc0b
using Base.Threads

# ╔═╡ 3c21884c-c58c-4400-a16b-8a6a97305a4e
pyplot()

# ╔═╡ fa810d48-c95f-43af-8378-ecd0558ff253
hk(n::Integer,x::Float64)=sf_bessel_jl(n,x)+im*sf_bessel_yl(n,x);

# ╔═╡ 348d131d-73c4-4476-a866-a511e36cc8bf
function ymn(n::Integer,m::Integer,θ::Float64,ϕ::Float64) 	#spherical harmonics
if m>=0
	return sf_legendre_sphPlm(n,m,cos(θ))*exp(im*m*ϕ);
elseif isodd(-m)
	return -sf_legendre_sphPlm(n,-m,cos(θ))*exp(im*m*ϕ);
else
	return sf_legendre_sphPlm(n,-m,cos(θ))*exp(im*m*ϕ);
end
end

# ╔═╡ 3ec16d5d-2d1d-4de3-864c-e7d452bd3af0
##########################build simulation model###############################
function buildModelProto(dimensionTmp::Integer, ρbTmp::Float64, cbTmp::Float64, ρpTmp::Float64, cpTmp::Float64, ωbTmp::Float64, pNoTmp::Integer, RsTmp::Float64, parPos::Matrix{Float64})
parPosTmp=deepcopy(parPos);
######################set incident waves#######################################
incAmpTmp = 3000.0;	#amplitude of incident beam
incDirTmp1 = [0.0, 0.0, 1.0];	#set incident direction 1
incDirTmp2 = [0.0, 0.0, -1.0];	#set incident direction 2
incPosTmp1 = [0.0, 0.0, -1];	#original position of incident wave 1
incPosTmp2 = [0.0, 0.0, 1];	#original position of incident wave 2
bgMediumTmp = Acoustic(dimensionTmp; ρ = ρbTmp, c = cbTmp);	#build background acoustic model
waveTmp = plane_source(bgMediumTmp; amplitude = incAmpTmp, direction = incDirTmp1, position = incPosTmp1)+plane_source(bgMediumTmp; amplitude = -incAmpTmp, direction = incDirTmp2, position = incPosTmp2);	#build incident plane wave
#########################incident wave done########################################
#########################set particles#############################################
parMediumTmp = Acoustic(dimensionTmp; ρ = ρpTmp, c = cpTmp);	#build the acoustic model in particles
particlesTmp=Array{Particle{dimensionTmp, Acoustic{Float64, dimensionTmp}, Sphere{Float64, dimensionTmp}}}(undef, 0);	#define a null array to store particles model
#build particle set
for iTmp in 1:pNoTmp
	parShapeTmp=Sphere(parPosTmp[iTmp,:],RsTmp);
	particlesTmp=push!(particlesTmp,Particle(parMediumTmp,parShapeTmp));
end
##########################particles done############################################
simModelTmp=FrequencySimulation(particlesTmp,waveTmp);#build simulation model
return simModelTmp
end

# ╔═╡ f8647bf0-9b6a-42b2-b55a-177c6d828628
########################get expansion coefficients###########################
function getCoefProto(ωbTmp::Float64, modelTmp::FrequencySimulation, coefOrderTmp::Integer)
simModelTmp=modelTmp;
coefDataTmp=basis_coefficients(simModelTmp,ωbTmp,basis_order=coefOrderTmp);#store the expansion coefficients
return coefDataTmp
end

# ╔═╡ 333ce538-8a65-4bcc-926f-fa99dd3692d4
###########define function to calculate pressure in position [x,y,z]############
function pProto(x::Float64, y::Float64, z::Float64, ωbTmp::Float64, modelTmp::FrequencySimulation, coefData::Matrix{ComplexF64})
coefDataTmp=deepcopy(coefData);
simModelTmp=modelTmp;
pNoTmp=length(simModelTmp.particles);
dimensionTmp=typeof(simModelTmp.source.medium).parameters[2];
cbTmp=Float64(simModelTmp.source.medium.c);     #make soundspeed a real number
coefOrderTmp=Int(sqrt(length(coefDataTmp[:,1])))-1;
parPosTmp=Matrix{Float64}(undef,pNoTmp,dimensionTmp);
k=ωbTmp/cbTmp;
for iTmp in 1:pNoTmp
	for jTmp in 1:dimensionTmp
		parPosTmp[iTmp,jTmp]=simModelTmp.particles[iTmp].shape.origin[jTmp];
	end
end
pField=0.0+0.0*im;
r=Array{Float64}(undef,pNoTmp);
θ=Array{Float64}(undef,pNoTmp);
ϕ=Array{Float64}(undef,pNoTmp);
for iTmp in 1:pNoTmp
	parPosTmp2=deepcopy(parPosTmp[iTmp,:]);
	xx=x-parPosTmp2[1];
	yy=y-parPosTmp2[2];
	zz=z-parPosTmp2[3];
	r[iTmp]=sqrt(xx*xx+yy*yy+zz*zz);
	θ[iTmp]=acos(zz/max(r[iTmp],0.000000001));
	if yy==0&&xx>0
		ϕ[iTmp]=0.0;
	elseif yy==0&&xx<0
		ϕ[iTmp]=π;
	else
		ϕ[iTmp]=(1-sign(yy))*π+sign(yy)*acos(xx/max(sqrt(xx*xx+yy*yy),0.00000001));
	end
	for nTmp in 0:coefOrderTmp
		for mTmp in -nTmp:nTmp
			pField+=coefDataTmp[nTmp*nTmp+nTmp+mTmp+1,iTmp]*hk(nTmp,k*r[iTmp])*ymn(nTmp,mTmp,θ[iTmp],ϕ[iTmp]);
		end
	end
end
return pField+simModelTmp.source.field([x,y,z],ωbTmp);
end

# ╔═╡ 685c0346-b182-4491-98c0-83457aa99528
############Calculate velocity field by using five-point stencil###############
function vProto(x::Float64, y::Float64, z::Float64, ωbTmp::Float64, modelTmp::FrequencySimulation, coefData::Matrix{ComplexF64})
coefDataTmp=deepcopy(coefData);
simModelTmp=modelTmp;
dimensionTmp=typeof(simModelTmp.source.medium).parameters[2];
δh=0.00005;
ρbTmp=simModelTmp.source.medium.ρ;
vField=Array{ComplexF64}(undef,dimensionTmp);
pX2=pProto(x+2*δh,y,z,ωbTmp,simModelTmp,coefDataTmp);
pX1=pProto(x+δh,y,z,ωbTmp,simModelTmp,coefDataTmp);
pXN2=pProto(x-2*δh,y,z,ωbTmp,simModelTmp,coefDataTmp);	#N denotes negative
pXN1=pProto(x-δh,y,z,ωbTmp,simModelTmp,coefDataTmp);
pY2=pProto(x,y+2*δh,z,ωbTmp,simModelTmp,coefDataTmp);
pY1=pProto(x,y+δh,z,ωbTmp,simModelTmp,coefDataTmp);
pYN2=pProto(x,y-2*δh,z,ωbTmp,simModelTmp,coefDataTmp);
pYN1=pProto(x,y-δh,z,ωbTmp,simModelTmp,coefDataTmp);
pZ2=pProto(x,y,z+2*δh,ωbTmp,simModelTmp,coefDataTmp);
pZ1=pProto(x,y,z+δh,ωbTmp,simModelTmp,coefDataTmp);
pZN2=pProto(x,y,z-2*δh,ωbTmp,simModelTmp,coefDataTmp);
pZN1=pProto(x,y,z-δh,ωbTmp,simModelTmp,coefDataTmp);
vField[1]=-im/ρbTmp/ωbTmp*(-pX2+8*pX1-8*pXN1+pXN2)/12/δh; 
vField[2]=-im/ρbTmp/ωbTmp*(-pY2+8*pY1-8*pYN1+pYN2)/12/δh; 
vField[3]=-im/ρbTmp/ωbTmp*(-pZ2+8*pZ1-8*pZN1+pZN2)/12/δh; 
return vField
end

# ╔═╡ 6230e3b7-1748-4f4b-a0b7-364389526809
############calculate time average of stress tensor#################
function T(x::Float64, y::Float64, z::Float64, ωbTmp::Float64, modelTmp::FrequencySimulation, coefData::Matrix{ComplexF64})
coefDataTmp=deepcopy(coefData);
simModelTmp=modelTmp;
dimensionTmp=typeof(simModelTmp.source.medium).parameters[2];
TDataTmp=Matrix{Float64}(undef,dimensionTmp,dimensionTmp);	#matrix to store stress tensor in the form:
##############################
######### Txx Txy Txz ########
######### Tyx Tyy Tyz ########
######### Tzx Tzy Tzz ########
##############################
ρbTmp=simModelTmp.source.medium.ρ;
cbTmp=simModelTmp.source.medium.c;
vTmp=vProto(x,y,z,ωbTmp,simModelTmp,coefDataTmp)
vxTmp=vTmp[1];
vyTmp=vTmp[2];
vzTmp=vTmp[3];
cjvxTmp=conj(vxTmp);
cjvyTmp=conj(vyTmp);
cjvzTmp=conj(vzTmp);
pTmp=pProto(x,y,z,ωbTmp,simModelTmp,coefDataTmp);
cjpTmp=conj(pTmp);
vSqu=real(vxTmp*cjvxTmp+vyTmp*cjvyTmp+vzTmp*cjvzTmp);
pSquCoeff=real(pTmp*cjpTmp/2/ρbTmp/cbTmp/cbTmp);
TDataTmp[1,1]=0.5*(ρbTmp*(real(vxTmp*cjvxTmp)-0.5*vSqu)+pSquCoeff);
TDataTmp[1,2]=0.5*ρbTmp*real(vxTmp*cjvyTmp);
TDataTmp[1,3]=0.5*ρbTmp*real(vxTmp*cjvzTmp);
TDataTmp[2,1]=0.5*ρbTmp*real(vyTmp*cjvxTmp);
TDataTmp[2,2]=0.5*(ρbTmp*(real(vyTmp*cjvyTmp)-0.5*vSqu)+pSquCoeff);
TDataTmp[2,3]=0.5*ρbTmp*real(vyTmp*cjvzTmp);
TDataTmp[3,1]=0.5*ρbTmp*real(vzTmp*cjvxTmp);
TDataTmp[3,2]=0.5*ρbTmp*real(vzTmp*cjvyTmp);
TDataTmp[3,3]=0.5*(ρbTmp*(real(vzTmp*cjvzTmp)-0.5*vSqu)+pSquCoeff);
return TDataTmp
end

# ╔═╡ 5615d75d-e805-4bc6-b5f2-4cdd278dd996
###################calculate force density##########################
function fDens(θ::Float64,ϕ::Float64,parID::Integer,ωbTmp::Float64, modelTmp::FrequencySimulation, coefData::Matrix{ComplexF64})
coefDataTmp=deepcopy(coefData);
simModelTmp=modelTmp;
dimensionTmp=typeof(simModelTmp.source.medium).parameters[2];
fDensDataTmp=Array{Float64}(undef,dimensionTmp);	#[fDensx,fDensy,fDensz]
R=simModelTmp.particles[parID].shape.radius+0.00015;
x0=simModelTmp.particles[parID].shape.origin[1];#x0,y0,z0 denote particle's position
y0=simModelTmp.particles[parID].shape.origin[2];
z0=simModelTmp.particles[parID].shape.origin[3];
sθ=sin(θ);cθ=cos(θ);sϕ=sin(ϕ);cϕ=cos(ϕ);
xTmp=x0+R*sθ*cϕ;
yTmp=y0+R*sθ*sϕ;
zTmp=z0+R*cθ;
TData=T(xTmp,yTmp,zTmp,ωbTmp,simModelTmp,coefDataTmp);
fDensDataTmp[1]=(TData[1,1]*sθ*cϕ+TData[1,2]*sθ*sϕ+TData[1,3]*cθ)*R*R*sθ;
fDensDataTmp[2]=(TData[2,1]*sθ*cϕ+TData[2,2]*sθ*sϕ+TData[2,3]*cθ)*R*R*sθ;
fDensDataTmp[3]=(TData[3,1]*sθ*cϕ+TData[3,2]*sθ*sϕ+TData[3,3]*cθ)*R*R*sθ;
return fDensDataTmp
end

# ╔═╡ e7dd7364-0b1c-45bb-b156-16a07ff5580e
###########function to execute part of the parallel computation###############
function forcePart1(parID::Integer,n1::Integer,n2::Integer,lowLim1::Float64,upLim1::Float64,lowLim2::Float64,upLim2::Float64,ωbTmp::Float64, modelTmp::FrequencySimulation, coefData::Matrix{ComplexF64})
coefDataTmp=deepcopy(coefData);
fDataTmp=[0.0,0.0,0.0]############[fx,fy,fz]#######################
xDataStoreTmp=zeros(n1,n2);
yDataStoreTmp=zeros(n1,n2);
zDataStoreTmp=zeros(n1,n2);
step1=(upLim1-lowLim1)/n1;
step2=(upLim2-lowLim2)/n2;
for i in 1:n1
	θTmp=lowLim1+(i-1)*step1;
	for j in 1:n2
		ϕTmp=lowLim2+(j-1)*step2;
		fDensData=fDens(θTmp,ϕTmp,parID,ωbTmp,modelTmp,coefDataTmp);
		xDataStoreTmp[i,j]=fDensData[1];
		yDataStoreTmp[i,j]=fDensData[2];
		zDataStoreTmp[i,j]=fDensData[3];
	end
end
fDataTmp[1]=-sum(xDataStoreTmp)*step1*step2;
fDataTmp[2]=-sum(yDataStoreTmp)*step1*step2;
fDataTmp[3]=-sum(zDataStoreTmp)*step1*step2;
return fDataTmp
end

# ╔═╡ 5eb23234-6a70-4aa7-b569-816b05e84dd9
###########function to execute part of the parallel computation###############
function forcePart2(parID::Integer,n1::Integer,n2::Integer,lowLim1::Float64,upLim1::Float64,lowLim2::Float64,upLim2::Float64,ωbTmp::Float64, modelTmp::FrequencySimulation, coefData::Matrix{ComplexF64})
coefDataTmp=deepcopy(coefData);
fDataTmp=[0.0,0.0,0.0]############[fx,fy,fz]#######################
xDataStoreTmp=zeros(n1,n2);
yDataStoreTmp=zeros(n1,n2);
zDataStoreTmp=zeros(n1,n2);
step1=(upLim1-lowLim1)/n1;
step2=(upLim2-lowLim2)/n2;
for i in 1:n1
	θTmp=lowLim1+(i-1)*step1;
	for j in 1:n2
		ϕTmp=lowLim2+(j-1)*step2;
		fDensData=fDens(θTmp,ϕTmp,parID,ωbTmp,modelTmp,coefDataTmp);
		xDataStoreTmp[i,j]=fDensData[1];
		yDataStoreTmp[i,j]=fDensData[2];
		zDataStoreTmp[i,j]=fDensData[3];
	end
end
fDataTmp[1]=-sum(xDataStoreTmp)*step1*step2;
fDataTmp[2]=-sum(yDataStoreTmp)*step1*step2;
fDataTmp[3]=-sum(zDataStoreTmp)*step1*step2;
return fDataTmp
end

# ╔═╡ 2b76f2d5-4194-4933-a834-52f9c1a4de5a
###########function to execute part of the parallel computation###############
function forcePart3(parID::Integer,n1::Integer,n2::Integer,lowLim1::Float64,upLim1::Float64,lowLim2::Float64,upLim2::Float64,ωbTmp::Float64, modelTmp::FrequencySimulation, coefData::Matrix{ComplexF64})
coefDataTmp=deepcopy(coefData);
fDataTmp=[0.0,0.0,0.0]############[fx,fy,fz]#######################
xDataStoreTmp=zeros(n1,n2);
yDataStoreTmp=zeros(n1,n2);
zDataStoreTmp=zeros(n1,n2);
step1=(upLim1-lowLim1)/n1;
step2=(upLim2-lowLim2)/n2;
for i in 1:n1
	θTmp=lowLim1+(i-1)*step1;
	for j in 1:n2
		ϕTmp=lowLim2+(j-1)*step2;
		fDensData=fDens(θTmp,ϕTmp,parID,ωbTmp,modelTmp,coefDataTmp);
		xDataStoreTmp[i,j]=fDensData[1];
		yDataStoreTmp[i,j]=fDensData[2];
		zDataStoreTmp[i,j]=fDensData[3];
	end
end
fDataTmp[1]=-sum(xDataStoreTmp)*step1*step2;
fDataTmp[2]=-sum(yDataStoreTmp)*step1*step2;
fDataTmp[3]=-sum(zDataStoreTmp)*step1*step2;
return fDataTmp
end

# ╔═╡ e6d5b885-146d-44d1-bc9e-78b6bcadebb0
###########function to execute part of the parallel computation###############
function forcePart4(parID::Integer,n1::Integer,n2::Integer,lowLim1::Float64,upLim1::Float64,lowLim2::Float64,upLim2::Float64,ωbTmp::Float64, modelTmp::FrequencySimulation, coefData::Matrix{ComplexF64})
coefDataTmp=deepcopy(coefData);
fDataTmp=[0.0,0.0,0.0]############[fx,fy,fz]#######################
xDataStoreTmp=zeros(n1,n2);
yDataStoreTmp=zeros(n1,n2);
zDataStoreTmp=zeros(n1,n2);
step1=(upLim1-lowLim1)/n1;
step2=(upLim2-lowLim2)/n2;
for i in 1:n1
	θTmp=lowLim1+(i-1)*step1;
	for j in 1:n2
		ϕTmp=lowLim2+(j-1)*step2;
		fDensData=fDens(θTmp,ϕTmp,parID,ωbTmp,modelTmp,coefDataTmp);
		xDataStoreTmp[i,j]=fDensData[1];
		yDataStoreTmp[i,j]=fDensData[2];
		zDataStoreTmp[i,j]=fDensData[3];
	end
end
fDataTmp[1]=-sum(xDataStoreTmp)*step1*step2;
fDataTmp[2]=-sum(yDataStoreTmp)*step1*step2;
fDataTmp[3]=-sum(zDataStoreTmp)*step1*step2;
return fDataTmp
end

# ╔═╡ 5e4b3511-0132-4c7a-987b-70e57a82a48f
###########function to execute part of the parallel computation###############
function forcePart5(parID::Integer,n1::Integer,n2::Integer,lowLim1::Float64,upLim1::Float64,lowLim2::Float64,upLim2::Float64,ωbTmp::Float64, modelTmp::FrequencySimulation, coefData::Matrix{ComplexF64})
coefDataTmp=deepcopy(coefData);
fDataTmp=[0.0,0.0,0.0]############[fx,fy,fz]#######################
xDataStoreTmp=zeros(n1,n2);
yDataStoreTmp=zeros(n1,n2);
zDataStoreTmp=zeros(n1,n2);
step1=(upLim1-lowLim1)/n1;
step2=(upLim2-lowLim2)/n2;
for i in 1:n1
	θTmp=lowLim1+(i-1)*step1;
	for j in 1:n2
		ϕTmp=lowLim2+(j-1)*step2;
		fDensData=fDens(θTmp,ϕTmp,parID,ωbTmp,modelTmp,coefDataTmp);
		xDataStoreTmp[i,j]=fDensData[1];
		yDataStoreTmp[i,j]=fDensData[2];
		zDataStoreTmp[i,j]=fDensData[3];
	end
end
fDataTmp[1]=-sum(xDataStoreTmp)*step1*step2;
fDataTmp[2]=-sum(yDataStoreTmp)*step1*step2;
fDataTmp[3]=-sum(zDataStoreTmp)*step1*step2;
return fDataTmp
end

# ╔═╡ ebe79ddc-137c-4070-83da-e1c91e179c90
###########function to execute part of the parallel computation###############
function forcePart6(parID::Integer,n1::Integer,n2::Integer,lowLim1::Float64,upLim1::Float64,lowLim2::Float64,upLim2::Float64,ωbTmp::Float64, modelTmp::FrequencySimulation, coefData::Matrix{ComplexF64})
coefDataTmp=deepcopy(coefData);
fDataTmp=[0.0,0.0,0.0]############[fx,fy,fz]#######################
xDataStoreTmp=zeros(n1,n2);
yDataStoreTmp=zeros(n1,n2);
zDataStoreTmp=zeros(n1,n2);
step1=(upLim1-lowLim1)/n1;
step2=(upLim2-lowLim2)/n2;
for i in 1:n1
	θTmp=lowLim1+(i-1)*step1;
	for j in 1:n2
		ϕTmp=lowLim2+(j-1)*step2;
		fDensData=fDens(θTmp,ϕTmp,parID,ωbTmp,modelTmp,coefDataTmp);
		xDataStoreTmp[i,j]=fDensData[1];
		yDataStoreTmp[i,j]=fDensData[2];
		zDataStoreTmp[i,j]=fDensData[3];
	end
end
fDataTmp[1]=-sum(xDataStoreTmp)*step1*step2;
fDataTmp[2]=-sum(yDataStoreTmp)*step1*step2;
fDataTmp[3]=-sum(zDataStoreTmp)*step1*step2;
return fDataTmp
end

# ╔═╡ b058329a-10ab-48b9-b51e-272982a02c17
###########function to execute part of the parallel computation###############
function forcePart7(parID::Integer,n1::Integer,n2::Integer,lowLim1::Float64,upLim1::Float64,lowLim2::Float64,upLim2::Float64,ωbTmp::Float64, modelTmp::FrequencySimulation, coefData::Matrix{ComplexF64})
coefDataTmp=deepcopy(coefData);
fDataTmp=[0.0,0.0,0.0]############[fx,fy,fz]#######################
xDataStoreTmp=zeros(n1,n2);
yDataStoreTmp=zeros(n1,n2);
zDataStoreTmp=zeros(n1,n2);
step1=(upLim1-lowLim1)/n1;
step2=(upLim2-lowLim2)/n2;
for i in 1:n1
	θTmp=lowLim1+(i-1)*step1;
	for j in 1:n2
		ϕTmp=lowLim2+(j-1)*step2;
		fDensData=fDens(θTmp,ϕTmp,parID,ωbTmp,modelTmp,coefDataTmp);
		xDataStoreTmp[i,j]=fDensData[1];
		yDataStoreTmp[i,j]=fDensData[2];
		zDataStoreTmp[i,j]=fDensData[3];
	end
end
fDataTmp[1]=-sum(xDataStoreTmp)*step1*step2;
fDataTmp[2]=-sum(yDataStoreTmp)*step1*step2;
fDataTmp[3]=-sum(zDataStoreTmp)*step1*step2;
return fDataTmp
end

# ╔═╡ fa6e0895-8de3-48ee-b1bf-39974d9a843a
###########function to execute part of the parallel computation###############
function forcePart8(parID::Integer,n1::Integer,n2::Integer,lowLim1::Float64,upLim1::Float64,lowLim2::Float64,upLim2::Float64,ωbTmp::Float64, modelTmp::FrequencySimulation, coefData::Matrix{ComplexF64})
coefDataTmp=deepcopy(coefData);
fDataTmp=[0.0,0.0,0.0]############[fx,fy,fz]#######################
xDataStoreTmp=zeros(n1,n2);
yDataStoreTmp=zeros(n1,n2);
zDataStoreTmp=zeros(n1,n2);
step1=(upLim1-lowLim1)/n1;
step2=(upLim2-lowLim2)/n2;
for i in 1:n1
	θTmp=lowLim1+(i-1)*step1;
	for j in 1:n2
		ϕTmp=lowLim2+(j-1)*step2;
		fDensData=fDens(θTmp,ϕTmp,parID,ωbTmp,modelTmp,coefDataTmp);
		xDataStoreTmp[i,j]=fDensData[1];
		yDataStoreTmp[i,j]=fDensData[2];
		zDataStoreTmp[i,j]=fDensData[3];
	end
end
fDataTmp[1]=-sum(xDataStoreTmp)*step1*step2;
fDataTmp[2]=-sum(yDataStoreTmp)*step1*step2;
fDataTmp[3]=-sum(zDataStoreTmp)*step1*step2;
return fDataTmp
end

# ╔═╡ 0844d606-d104-44ce-bf09-967fbf708629
##################single-threaded routine#################################
################numerical quadrature for acoustic force################
#############n1,n2 denote sample numbers in polar angle θ and azimuthal angle ϕ coordinate respectively###############
function forceT1(parID::Integer,n1::Integer,n2::Integer,ωbTmp::Float64, modelTmp::FrequencySimulation, coefData::Matrix{ComplexF64})
coefDataTmp=deepcopy(coefData);
fDataTmp=[0.0,0.0,0.0]############[fx,fy,fz]#######################
xDataStoreTmp=zeros(n1,n2);
yDataStoreTmp=zeros(n1,n2);
zDataStoreTmp=zeros(n1,n2);
step1=π/n1;
step2=2*π/n2;
for i in 1:n1
	θTmp=(i-1)*step1;
	for j in 1:n2
		ϕTmp=(j-1)*step2;
		fDensData=fDens(θTmp,ϕTmp,parID,ωbTmp,modelTmp,coefDataTmp);
		xDataStoreTmp[i,j]=fDensData[1];
		yDataStoreTmp[i,j]=fDensData[2];
		zDataStoreTmp[i,j]=fDensData[3];
	end
end
fDataTmp[1]=-sum(xDataStoreTmp)*step1*step2;
fDataTmp[2]=-sum(yDataStoreTmp)*step1*step2;
fDataTmp[3]=-sum(zDataStoreTmp)*step1*step2;
return fDataTmp
end

# ╔═╡ e5033c56-20bd-42e1-bd24-4dda795d5f83
##################2-threaded routine#################################
################numerical quadrature for acoustic force################
#############n1,n2 denote sample numbers in polar angle θ and azimuthal angle ϕ coordinate respectively###############
function forceT2(parID::Integer,n1::Integer,n2::Integer,ωbTmp::Float64, modelTmp::FrequencySimulation, coefData::Matrix{ComplexF64})
threadNo=2;
threadNo1=1;
threadNo2=2;
n1Tmp=Int(n1/threadNo1);
n2Tmp=Int(n2/threadNo2);
lowLim1=Array{Float64}(undef,threadNo1);
upLim1=Array{Float64}(undef,threadNo1);
lowLim2=Array{Float64}(undef,threadNo2);
upLim2=Array{Float64}(undef,threadNo2);
for i in 1:threadNo1
	lowLim1[i]=(i-1)*π/threadNo1;
	upLim1[i]=lowLim1[i]+π/threadNo1;
end
for i in 1:threadNo2
	lowLim2[i]=(i-1)*2π/threadNo2;
	upLim2[i]=lowLim2[i]+2π/threadNo2;
end
result1=Threads.@spawn 	forcePart1(parID,n1Tmp,n2Tmp,lowLim1[1],upLim1[1],lowLim2[1],upLim2[1],ωbTmp, modelTmp, coefData);
result2=Threads.@spawn 	forcePart2(parID,n1Tmp,n2Tmp,lowLim1[1],upLim1[1],lowLim2[2],upLim2[2],ωbTmp, modelTmp, coefData);
return fetch(result1) + fetch(result2)
end

# ╔═╡ 6556540d-19bb-405f-a191-8bbcd63e224e
##################4-threaded routine#################################
################numerical quadrature for acoustic force################
#############n1,n2 denote sample numbers in polar angle θ and azimuthal angle ϕ coordinate respectively###############
function forceT4(parID::Integer,n1::Integer,n2::Integer,ωbTmp::Float64, modelTmp::FrequencySimulation, coefData::Matrix{ComplexF64})
threadNo=4;
threadNo1=2;
threadNo2=2;
n1Tmp=Int(n1/threadNo1);
n2Tmp=Int(n2/threadNo2);
lowLim1=Array{Float64}(undef,threadNo1);
upLim1=Array{Float64}(undef,threadNo1);
lowLim2=Array{Float64}(undef,threadNo2);
upLim2=Array{Float64}(undef,threadNo2);
for i in 1:threadNo1
	lowLim1[i]=(i-1)*π/threadNo1;
	upLim1[i]=lowLim1[i]+π/threadNo1;
end
for i in 1:threadNo2
	lowLim2[i]=(i-1)*2π/threadNo2;
	upLim2[i]=lowLim2[i]+2π/threadNo2;
end
result1=Threads.@spawn 	forcePart1(parID,n1Tmp,n2Tmp,lowLim1[1],upLim1[1],lowLim2[1],upLim2[1],ωbTmp, modelTmp, coefData);
result2=Threads.@spawn 	forcePart2(parID,n1Tmp,n2Tmp,lowLim1[1],upLim1[1],lowLim2[2],upLim2[2],ωbTmp, modelTmp, coefData);
result3=Threads.@spawn 	forcePart3(parID,n1Tmp,n2Tmp,lowLim1[2],upLim1[2],lowLim2[1],upLim2[1],ωbTmp, modelTmp, coefData);
result4=Threads.@spawn 	forcePart4(parID,n1Tmp,n2Tmp,lowLim1[2],upLim1[2],lowLim2[2],upLim2[2],ωbTmp, modelTmp, coefData);
return fetch(result1) + fetch(result2) + fetch(result3) + fetch(result4)
end

# ╔═╡ 48600999-3c69-48fb-b827-90d7e655e18f
##################6-threaded routine#################################
################numerical quadrature for acoustic force################
#############n1,n2 denote sample numbers in polar angle θ and azimuthal angle ϕ coordinate respectively###############
function forceT6(parID::Integer,n1::Integer,n2::Integer,ωbTmp::Float64, modelTmp::FrequencySimulation, coefData::Matrix{ComplexF64})
threadNo=4;
threadNo1=2;
threadNo2=3;
n1Tmp=Int(n1/threadNo1);
n2Tmp=Int(n2/threadNo2);
lowLim1=Array{Float64}(undef,threadNo1);
upLim1=Array{Float64}(undef,threadNo1);
lowLim2=Array{Float64}(undef,threadNo2);
upLim2=Array{Float64}(undef,threadNo2);
for i in 1:threadNo1
	lowLim1[i]=(i-1)*π/threadNo1;
	upLim1[i]=lowLim1[i]+π/threadNo1;
end
for i in 1:threadNo2
	lowLim2[i]=(i-1)*2π/threadNo2;
	upLim2[i]=lowLim2[i]+2π/threadNo2;
end
result1=Threads.@spawn 	forcePart1(parID,n1Tmp,n2Tmp,lowLim1[1],upLim1[1],lowLim2[1],upLim2[1],ωbTmp, modelTmp, coefData);
result2=Threads.@spawn 	forcePart2(parID,n1Tmp,n2Tmp,lowLim1[1],upLim1[1],lowLim2[2],upLim2[2],ωbTmp, modelTmp, coefData);
result3=Threads.@spawn 	forcePart3(parID,n1Tmp,n2Tmp,lowLim1[1],upLim1[1],lowLim2[3],upLim2[3],ωbTmp, modelTmp, coefData);
result4=Threads.@spawn 	forcePart4(parID,n1Tmp,n2Tmp,lowLim1[2],upLim1[2],lowLim2[1],upLim2[1],ωbTmp, modelTmp, coefData);
result5=Threads.@spawn 	forcePart5(parID,n1Tmp,n2Tmp,lowLim1[2],upLim1[2],lowLim2[2],upLim2[2],ωbTmp, modelTmp, coefData);
result6=Threads.@spawn 	forcePart6(parID,n1Tmp,n2Tmp,lowLim1[2],upLim1[2],lowLim2[3],upLim2[3],ωbTmp, modelTmp, coefData);
return fetch(result1) + fetch(result2) + fetch(result3) + fetch(result4)+ fetch(result5) + fetch(result6)
end

# ╔═╡ 82f0e9b0-6268-4b1c-99f2-402e2a25cdf3
##################8-threaded routine#################################
################numerical quadrature for acoustic force################
#############n1,n2 denote sample numbers in polar angle θ and azimuthal angle ϕ coordinate respectively###############
function forceT8(parID::Integer,n1::Integer,n2::Integer,ωbTmp::Float64, modelTmp::FrequencySimulation, coefData::Matrix{ComplexF64})
threadNo=4;
threadNo1=2;
threadNo2=4;
n1Tmp=Int(n1/threadNo1);
n2Tmp=Int(n2/threadNo2);
lowLim1=Array{Float64}(undef,threadNo1);
upLim1=Array{Float64}(undef,threadNo1);
lowLim2=Array{Float64}(undef,threadNo2);
upLim2=Array{Float64}(undef,threadNo2);
for i in 1:threadNo1
	lowLim1[i]=(i-1)*π/threadNo1;
	upLim1[i]=lowLim1[i]+π/threadNo1;
end
for i in 1:threadNo2
	lowLim2[i]=(i-1)*2π/threadNo2;
	upLim2[i]=lowLim2[i]+2π/threadNo2;
end
result1=Threads.@spawn 	forcePart1(parID,n1Tmp,n2Tmp,lowLim1[1],upLim1[1],lowLim2[1],upLim2[1],ωbTmp, modelTmp, coefData);
result2=Threads.@spawn 	forcePart2(parID,n1Tmp,n2Tmp,lowLim1[1],upLim1[1],lowLim2[2],upLim2[2],ωbTmp, modelTmp, coefData);
result3=Threads.@spawn 	forcePart3(parID,n1Tmp,n2Tmp,lowLim1[1],upLim1[1],lowLim2[3],upLim2[3],ωbTmp, modelTmp, coefData);
result4=Threads.@spawn 	forcePart4(parID,n1Tmp,n2Tmp,lowLim1[1],upLim1[1],lowLim2[4],upLim2[4],ωbTmp, modelTmp, coefData);
result5=Threads.@spawn 	forcePart5(parID,n1Tmp,n2Tmp,lowLim1[2],upLim1[2],lowLim2[1],upLim2[1],ωbTmp, modelTmp, coefData);
result6=Threads.@spawn 	forcePart6(parID,n1Tmp,n2Tmp,lowLim1[2],upLim1[2],lowLim2[2],upLim2[2],ωbTmp, modelTmp, coefData);
result7=Threads.@spawn 	forcePart7(parID,n1Tmp,n2Tmp,lowLim1[2],upLim1[2],lowLim2[3],upLim2[3],ωbTmp, modelTmp, coefData);
result8=Threads.@spawn 	forcePart8(parID,n1Tmp,n2Tmp,lowLim1[2],upLim1[2],lowLim2[4],upLim2[4],ωbTmp, modelTmp, coefData);
return fetch(result1) + fetch(result2) + fetch(result3) + fetch(result4)+ fetch(result5) + fetch(result6)+ fetch(result7) + fetch(result8)
end

# ╔═╡ 9d15bca5-8c13-4998-83b6-6cf49bf9e2f4
function force(parID::Integer,n1::Integer,n2::Integer,ωbTmp::Float64, modelTmp::FrequencySimulation, coefData::Matrix{ComplexF64})
coefDataTmp=deepcopy(coefData);
num_threads = Threads.nthreads();
if num_threads==1
	println("Single-threaded computation")
	return forceT1(parID,n1,n2,ωbTmp, modelTmp, coefDataTmp)
elseif 2<=num_threads&&num_threads<4
	println("2-threaded parallel computation")
	return forceT2(parID,n1,n2,ωbTmp, modelTmp, coefDataTmp)
elseif 4<=num_threads&&num_threads<6
	println("4-threaded parallel computation")
	return forceT4(parID,n1,n2,ωbTmp, modelTmp, coefDataTmp)
elseif 6<=num_threads&&num_threads<8
	println("6-threaded parallel computation")
	return forceT6(parID,n1,n2,ωbTmp, modelTmp, coefDataTmp)
else
	println("8-threaded parallel computation")
	return forceT8(parID,n1,n2,ωbTmp, modelTmp, coefDataTmp)
end
end


# ╔═╡ ff523df4-856c-42c5-ba13-9caf45ad8f32
###################calculate all particles forces###################
function allForce(n1::Integer,n2::Integer,ωbTmp::Float64, modelTmp::FrequencySimulation, coefData::Matrix{ComplexF64})
coefDataTmp=deepcopy(coefData);
simModelTmp=modelTmp;
pNoTmp=length(simModelTmp.particles);
dimensionTmp=typeof(simModelTmp.source.medium).parameters[2];
allForceTmp=Array{Float64}(undef,0)
for iTmp in 1:pNoTmp
	forceTmp=force(iTmp,n1,n2,ωbTmp, modelTmp, coefDataTmp);
	allForceTmp=push!(allForceTmp,forceTmp[1],forceTmp[2],forceTmp[3]);
end
return allForceTmp
end

# ╔═╡ e72bc97e-5c62-4fc5-b6ee-17651a37bedd
function forcePackLow(RsTmp::Float64,parPos::Matrix{Float64})
parPosTmp=deepcopy(parPos);
freqIn=40000;
ω=2.0*π*freqIn;
dimension=3;
ρb=1.225;
ρp=29.0;
cb=343.0;
cp=900.0;
pNo=length(parPosTmp[:,1]);
coefOrder=6;
modelTmp=buildModelProto(dimension, ρb, cb, ρp, cp, ω, pNo, RsTmp, parPosTmp)
coefData=getCoefProto(ω, modelTmp, coefOrder);
forceTmp=allForce(24,48,ω, modelTmp, coefData);
return forceTmp
end

# ╔═╡ 96de7809-1659-4ea4-bc53-9232d4e12270
function forcePackMiddle(RsTmp::Float64,parPos::Matrix{Float64})
parPosTmp=deepcopy(parPos);
freqIn=40000;
ω=2.0*π*freqIn;
dimension=3;
ρb=1.225;
ρp=29.0;
cb=343.0;
cp=900.0;
pNo=length(parPosTmp[:,1]);
coefOrder=8;
modelTmp=buildModelProto(dimension, ρb, cb, ρp, cp, ω, pNo, RsTmp, parPosTmp)
coefData=getCoefProto(ω, modelTmp, coefOrder);
forceTmp=allForce(30,60,ω, modelTmp, coefData);
return forceTmp
end

# ╔═╡ f386b9d4-50bf-4271-a167-953895164a32
function forcePackHigh(RsTmp::Float64,parPos::Matrix{Float64})
parPosTmp=deepcopy(parPos);
freqIn=40000;
ω=2.0*π*freqIn;
dimension=3;
ρb=1.225;
ρp=29.0;
cb=343.0;
cp=900.0;
pNo=length(parPosTmp[:,1]);
coefOrder=10;
modelTmp=buildModelProto(dimension, ρb, cb, ρp, cp, ω, pNo, RsTmp, parPosTmp)
coefData=getCoefProto(ω, modelTmp, coefOrder);
forceTmp=allForce(36,72,ω, modelTmp, coefData);
return forceTmp
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
GSL = "92c85e6c-cbff-5e0c-80f7-495c94daaecd"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
MultipleScattering = "80a8ab25-5750-5d93-a6d7-4adc97cdd5fb"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee"

[compat]
DelimitedFiles = "~1.9.1"
GSL = "~1.0.1"
MultipleScattering = "~0.1.16"
Plots = "~1.39.0"
PyPlot = "~2.11.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.2"
manifest_format = "2.0"
project_hash = "a360305cd01a43e938fc2f15e8dadbcdc56581c3"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "50171ec101d7e7557b1025425109727c84980a88"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.0"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "02aa26a4cf76381be7f66e020a3eddeb27b0a092"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "d9a8f86737b665e15a9641ecbac64deef9ce6724"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.23.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "e460f044ca8b99be31d35fe54fc33a5c33dd8ed7"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.9.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "5372dbbf8f0bdb8c700db5367132925c0771ef7e"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.2.1"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "8c86e48c0db1564a1d49548d3515ced5d604c408"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.9.1"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "e90caa41f5a86296e014e148ee061bd6c3edec96"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "d73afa4a2bb9de56077242d98cf763074ab9a970"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.9"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "1596bab77f4f073a14c62424283e7ebff3072eca"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.9+1"

[[deps.GSL]]
deps = ["GSL_jll", "Libdl", "Markdown"]
git-tree-sha1 = "3ebd07d519f5ec318d5bc1b4971e2472e14bd1f0"
uuid = "92c85e6c-cbff-5e0c-80f7-495c94daaecd"
version = "1.0.1"

[[deps.GSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "56f1e2c9e083e0bb7cf9a7055c280beb08a924c0"
uuid = "1b77fbbe-d8ee-58f0-85f9-836ddc23a7a4"
version = "2.7.2+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "cb56ccdd481c0dd7f975ad2b3b62d9eda088f7e2"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.9.14"

[[deps.HalfIntegers]]
git-tree-sha1 = "1cfb497b72e1e8ab2256334dee1aaad43fa279ad"
uuid = "f0d1745a-41c9-11e9-1dd9-e5d34d218721"
version = "1.5.1"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f689897ccbe049adb19a065c495e75f372ecd42b"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.4+0"

[[deps.LRUCache]]
git-tree-sha1 = "48c10e3cc27e30de82463c27bef0b8bdbd1dc634"
uuid = "8ac3fa9e-de4c-5943-b1dc-09c6b5f20637"
version = "1.4.1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "0d097476b6c381ab7906460ef1ef1638fbce1d91"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.2"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.MultipleScattering]]
deps = ["GSL", "LinearAlgebra", "OffsetArrays", "Printf", "ProgressMeter", "Random", "RecipesBase", "SpecialFunctions", "StaticArrays", "Statistics", "WignerSymbols"]
git-tree-sha1 = "cbab9acc6effc9c087202194a4b43e12e5dac9a9"
uuid = "80a8ab25-5750-5d93-a6d7-4adc97cdd5fb"
version = "0.1.16"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "2ac17d29c523ce1cd38e27785a7d23024853a4bb"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.10"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "bbb5c2115d63c2f1451cb70e5ef75e8fe4707019"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.22+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "ccee59c6e48e6f2edf8a5b64dc817b6729f99eb5"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.39.0"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "4c9f306e5d6603ae203c2000dd460d81a5251489"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.4"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "ae36206463b2395804f2787ffe172f44452b538d"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.8.0"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "43d304ac6f0354755f1d60730ece8c499980f7ba"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.96.1"

[[deps.PyPlot]]
deps = ["Colors", "LaTeXStrings", "PyCall", "Sockets", "Test", "VersionParsing"]
git-tree-sha1 = "9220a9dae0369f431168d60adab635f88aca7857"
uuid = "d330b81b-6aea-500a-939a-2ce795aea3ee"
version = "2.11.2"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "364898e8f13f7eaaceec55fd3d08680498c0aa6e"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.4.2+3"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RationalRoots]]
git-tree-sha1 = "52315cf3098691c1416a356925027af5ab5bf548"
uuid = "308eb6b3-cc68-5ff3-9e97-c3c4da4fa681"
version = "0.2.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "c60ec5c62180f27efea3ba2908480f8055e17cee"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "b0e70de949b55b1f8d2b9b69dcb378a02e5e9e00"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "0.12.6"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "75ebe04c5bed70b91614d684259b661c9e6274a4"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[deps.URIs]]
git-tree-sha1 = "b7a5e99f24892b6824a954199a45e9ffcc1c70f0"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "a72d22c7e13fe2de562feda8645aa134712a87ee"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.17.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "ed8d92d9774b077c53e1da50fd81a36af3744c1c"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.WignerSymbols]]
deps = ["HalfIntegers", "LRUCache", "Primes", "RationalRoots"]
git-tree-sha1 = "960e5f708871c1d9a28a7f1dbcaf4e0ee34ee960"
uuid = "9f57e263-0b3d-5e2e-b1be-24f2bb48858b"
version = "2.0.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "04a51d15436a572301b5abbb9d099713327e9fc4"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.4+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cf2c7de82431ca6f39250d2fc4aacd0daa1675c0"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╠═41b86c42-65a7-11ee-0959-ddd4537ec9be
# ╠═24689bb2-8862-442f-9b82-03b2270c3ed6
# ╠═4e14a103-6c4f-4bb0-8acc-c0be754be0c7
# ╠═57493ea3-4975-4ddf-8955-06ecbb1a03ff
# ╠═3c21884c-c58c-4400-a16b-8a6a97305a4e
# ╠═5bb8c76a-556b-4f2d-9af5-4c806220cba5
# ╠═ea5429e2-5e91-4003-bc4e-5ee1f1980056
# ╠═e245f8bd-f235-4ec6-b4e2-e405d650bc0b
# ╠═fa810d48-c95f-43af-8378-ecd0558ff253
# ╠═348d131d-73c4-4476-a866-a511e36cc8bf
# ╠═3ec16d5d-2d1d-4de3-864c-e7d452bd3af0
# ╠═f8647bf0-9b6a-42b2-b55a-177c6d828628
# ╠═333ce538-8a65-4bcc-926f-fa99dd3692d4
# ╠═685c0346-b182-4491-98c0-83457aa99528
# ╠═6230e3b7-1748-4f4b-a0b7-364389526809
# ╠═5615d75d-e805-4bc6-b5f2-4cdd278dd996
# ╠═e7dd7364-0b1c-45bb-b156-16a07ff5580e
# ╠═5eb23234-6a70-4aa7-b569-816b05e84dd9
# ╠═2b76f2d5-4194-4933-a834-52f9c1a4de5a
# ╠═e6d5b885-146d-44d1-bc9e-78b6bcadebb0
# ╠═5e4b3511-0132-4c7a-987b-70e57a82a48f
# ╠═ebe79ddc-137c-4070-83da-e1c91e179c90
# ╠═b058329a-10ab-48b9-b51e-272982a02c17
# ╠═fa6e0895-8de3-48ee-b1bf-39974d9a843a
# ╠═0844d606-d104-44ce-bf09-967fbf708629
# ╠═e5033c56-20bd-42e1-bd24-4dda795d5f83
# ╠═6556540d-19bb-405f-a191-8bbcd63e224e
# ╠═48600999-3c69-48fb-b827-90d7e655e18f
# ╠═82f0e9b0-6268-4b1c-99f2-402e2a25cdf3
# ╠═9d15bca5-8c13-4998-83b6-6cf49bf9e2f4
# ╠═ff523df4-856c-42c5-ba13-9caf45ad8f32
# ╠═e72bc97e-5c62-4fc5-b6ee-17651a37bedd
# ╠═96de7809-1659-4ea4-bc53-9232d4e12270
# ╠═f386b9d4-50bf-4271-a167-953895164a32
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002