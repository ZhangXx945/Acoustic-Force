### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 6f405d58-bdb6-48ff-8c87-b8eee3a854f5
#############generate particle position########################
function parPosGen(center::Bool,pNo::Integer,foldNo::Integer,ringNo::Integer,Rc::Array{Float64})
rcTmp=deepcopy(Rc);
if center
	parPosTmp[1,:]=[0.0 0 0];
	for i in 1:ringNo
		for j in 1:foldNo
			parPosTmp[j+(i-1)*foldNo+1,1]=rcTmp[i]*cos((j-1)*2π/foldNo);
			parPosTmp[j+(i-1)*foldNo+1,2]=rcTmp[i]*sin((j-1)*2π/foldNo);
			parPosTmp[j+(i-1)*foldNo+1,3]=0.0;
		end
	end
else
	for i in 1:ringNo
		for j in 1:foldNo
			parPosTmp[j+(i-1)*foldNo,1]=rcTmp[i]*cos((j-1)*2π/foldNo);
			parPosTmp[j+(i-1)*foldNo,2]=rcTmp[i]*sin((j-1)*2π/foldNo);
			parPosTmp[j+(i-1)*foldNo,3]=0.0;
		end
	end
end
	return parPosTmp
end

# ╔═╡ c930db94-1bc1-425b-9188-a1020a6b0b8c
##########molecule dynamics for equilibrium with rotational symmetries#############
function mdS(center::Bool,pNo::Integer,foldNo::Integer,ringNo::Integer,Rc::Array{Float64},Rs::Float64)
precSet=10^-6;	#set the precision of equilibrium position
rndDigNo=9;	#digit No. of round force.
minStepTmp=0.0002;	#minimal step
maxStepTmp=minStepTmp*10;	#maximal step
maxLopNo=500;	#max loop No to find equilibrium position
dimension=3;      #for 3D
ρp=29.0;         #particle density
parPosTmp=Matrix{Float64}(undef,pNo,3);
parPosData=Array{Float64}(undef,maxLopNo*10,dimension,pNo);
rcTmp=deepcopy(Rc);
parPosTmp=parPosGen(center,pNo,foldNo,ringNo,rcTmp);
forceTmp=round.(forcePackLow(Rs,parPosTmp);digits=rndDigNo);
disTmp2=round.(movDis(forceTmp,Rs,ρp);digits=rndDigNo);
if findmax(abs,disTmp2)[1]==0
	println("particles are in EQUILIBRIUM positions!")
	return parPosTmp
end
countTmp=1;
countTmp2=1;
while minStepTmp>=precSet
	minStepTmp=minStepTmp/5;
	maxStepTmp=maxStepTmp/5;
	for lopTmp in 1:maxLopNo
		disTmp=round.(movDis(forceTmp,Rs,ρp);digits=rndDigNo);
		if findmax(abs,disTmp)[1]==0
			println("particles are in EQUILIBRIUM positions!")
			return parPosTmp
		end
		if ensMov(disTmp)==1
			println("EQUILIBRIUM, ENSEMBLE MOVEMENT")
			println("Rs: "*string(Rs))
			println(parPosTmp)
			return parPosTmp
		end
		disTmp=disUpCheck(disTmp,maxStepTmp);
		disTmp=disLowCheck(disTmp,minStepTmp);
		disTmp=reScal(disTmp);
		if findmax(disTmp.*disTmp2)[1]<=0	#condition of equilibrium
			println("particles are in EQUILIBRIUM positions!")
			println("precision:",maxStepTmp)
			for i in 1:ringNo
				if center
					rcTmp[i]=rcTmp[i]+disTmp[(1+(i-1)*foldNo)*3+1]/2;
				else
					rcTmp[i]=rcTmp[i]+disTmp[(i-1)*foldNo*3+1]/2;
				end
			end
			parPosTmp=round.(parPosGen(center,pNo,foldNo,ringNo,rcTmp);digits=rndDigNo);
			disTmp2=deepcopy(disTmp);
			break
		else
			for i in 1:ringNo
				if center
					rcTmp[i]=rcTmp[i]+disTmp[(1+(i-1)*foldNo)*3+1,1];
				else
					rcTmp[i]=rcTmp[i]+disTmp[(i-1)*foldNo*3+1,1];
				end
			end
			parPosTmp=round.(parPosGen(center,pNo,foldNo,ringNo,rcTmp);digits=rndDigNo);
			for iTmp in 1:pNo
				for jTmp in 1:dimension
					parPosData[countTmp,jTmp,iTmp]=parPosTmp[iTmp,jTmp];
				end
			end
			disTmp2=deepcopy(disTmp);
			parPosTmp=sepCheck(parPosTmp,Rs);
			forceTmp=round.(forcePackLow(Rs,parPosTmp);digits=rndDigNo);
		end
		parPosData2=parPosData[1:countTmp,:,:];
		plot3d(plottitles=countTmp);
		println(["Max Displacement: " * string(maxStepTmp), "Total steps: " * string(countTmp),"RND No.: " * string(countTmp2),"Step No.: " * string(lopTmp)])
		for iTmp in 1:pNo
		display(path3d!(parPosData2[:,1,iTmp],parPosData2[:,2,iTmp],parPosData2[:,3,iTmp]))
		end
		println(disTmp)
		println("current particle positions:")
		println(parPosTmp)
		countTmp+=1;
	if lopTmp==maxLopNo
		println("Loop No. reaches Maximum!")
	end
	end
	countTmp2+=1;
end
return rcTmp
end


# ╔═╡ Cell order:
# ╠═6f405d58-bdb6-48ff-8c87-b8eee3a854f5
# ╠═c930db94-1bc1-425b-9188-a1020a6b0b8c