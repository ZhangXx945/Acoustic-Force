### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 6ae063fa-7b51-11ee-1653-d92b9fefdcdc
#################calculate force matrix######################
function forceMat(parPosTmp::Matrix{Float64},Rs::Float64)
parPos=deepcopy(parPosTmp);
δh=0.0001;
dimension=3;
dimension2=2;
rndDigNo=12;
pNo=length(parPos[:,1]);
kMat=Matrix{Float64}(undef,pNo*dimension2,pNo*dimension2);
for iTmp in 1:pNo
for jTmp in 1:dimension2
parPos[iTmp,jTmp]=parPos[iTmp,jTmp]+2*δh;
f1=round.(forcePackHigh(Rs,parPos); digits=rndDigNo);
parPos[iTmp,jTmp]=parPos[iTmp,jTmp]-2*δh;
parPos[iTmp,jTmp]=parPos[iTmp,jTmp]+δh;
f2=round.(forcePackHigh(Rs,parPos); digits=rndDigNo);
parPos[iTmp,jTmp]=parPos[iTmp,jTmp]-δh;
parPos[iTmp,jTmp]=parPos[iTmp,jTmp]-2*δh;
f3=round.(forcePackHigh(Rs,parPos); digits=rndDigNo);
parPos[iTmp,jTmp]=parPos[iTmp,jTmp]+2*δh;
parPos[iTmp,jTmp]=parPos[iTmp,jTmp]-δh;
f4=round.(forcePackHigh(Rs,parPos); digits=rndDigNo);
parPos[iTmp,jTmp]=parPos[iTmp,jTmp]+δh;
kk1=(-f1.+8*f2.-8*f4.+f3)/12/δh;
if dimension2==dimension
kk2=copy(kk1);
else
kk2=Array{Float64}(undef,pNo*dimension2);
countTmp=1;
for kTmp in 1:pNo*dimension
if kTmp%3==0
continue
else
kk2[countTmp]=kk1[kTmp];
countTmp+=1;
end
end
end
kMat[:,(iTmp-1)*dimension2+jTmp]=kk2;
end
println("particle "*string(iTmp)*" done")
end
return kMat
end

# ╔═╡ Cell order:
# ╠═6ae063fa-7b51-11ee-1653-d92b9fefdcdc