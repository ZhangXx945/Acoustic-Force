### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 6f405d58-bdb6-48ff-8c87-b8eee3a854f5
using MultipleScattering

# ╔═╡ ed3254e0-7b23-11ee-3e45-7dfc3a3015bf
##########generate movement according to force and particle property############
function movDis(fTmp2::Vector{Float64},Rs::Float64,ρp::Float64)
fTmp=deepcopy(fTmp2);
Δt=1;
aa=0.5*Δt*Δt*3/4/π/Rs/Rs/Rs/ρp;	#coefficient in calculating displacement in each step
return aa.*fTmp;
end

# ╔═╡ 05add776-b852-4e6f-9c2f-3ba7726617b0
###########################in case of overlapping##########################
function sepCheck(parPos::Matrix{Float64},Rs::Float64)
parPosTmp=deepcopy(parPos);
RTmp=Rs+0.00015;
pNo=length(parPosTmp[:,1]);
if pNo==1
	return parPosTmp
end
sepDisDataTmp=fill(16*RTmp*RTmp,sum(1:pNo-1));
countTmp=1;
for iTmp in 1:pNo-1
	for jTmp in iTmp+1:pNo
		sepDisDataTmp[countTmp]=sum((parPosTmp[iTmp,:].-parPosTmp[jTmp,:]).*(parPosTmp[iTmp,:].-parPosTmp[jTmp,:]));
		countTmp+=1;
	end
end
if findmin(sepDisDataTmp)[1]<=(2*RTmp)^2
	println("OVERLAPPING!")
	countTmp=1;
	for iTmp in 1:pNo-1
		for jTmp in iTmp+1:pNo
			if sepDisDataTmp[countTmp]<(2*RTmp)^2
				θTmp=acos((parPosTmp[jTmp,3]-parPosTmp[iTmp,3])/max(0.0000001,sqrt(sepDisDataTmp[countTmp])));
				yyTmp=parPosTmp[jTmp,2]-parPosTmp[iTmp,2];
				xxTmp=parPosTmp[jTmp,1]-parPosTmp[iTmp,1];
				if yyTmp==0&&xxTmp>0
					ϕTmp=0.0;
				elseif yyTmp==0&&xxTmp<0
					ϕTmp=π;
				else
					ϕTmp=(1-sign(yyTmp))*π+sign(yyTmp)*acos(xxTmp/max(sqrt(xxTmp*xxTmp+yyTmp*yyTmp),0.00000001));
				end
				parPosTmp[jTmp,1]=parPosTmp[iTmp,1]+3*Rs*sin(θTmp)*cos(ϕTmp);
				parPosTmp[jTmp,2]=parPosTmp[iTmp,2]+3*Rs*sin(θTmp)*sin(ϕTmp);
				#parPosTmp[jTmp,3]=parPosTmp[iTmp,3]+3*Rs*cos(θTmp);
				parPosTmp[jTmp,3]=0.0;
			end
			countTmp+=1;
		end
	end
	sepCheck(parPosTmp,Rs);
else
	return parPosTmp;
end
end

# ╔═╡ 64ddae5f-bea6-4ade-b4b7-ae66f43873f0
############transform displacement component into displacement###############
function disTrans(disTmp2::Vector{Float64})
	disTmp=deepcopy(disTmp2);
	dimension=3;     #for 3D
	pNo=Int(length(disTmp)/dimension);     #for 3D
	dis=Vector{Float64}(undef,pNo);
	for i in 1:pNo
		dis[i]=sqrt(disTmp[1+(i-1)*dimension]*disTmp[1+(i-1)*dimension]+disTmp[2+(i-1)*dimension]*disTmp[2+(i-1)*dimension]+disTmp[3+(i-1)*dimension]*disTmp[3+(i-1)*dimension]);
	end
	return dis
end

# ╔═╡ 9c25a47b-e3de-4747-b22a-40a10ded4668
#####################in case of too large movement#######################
function disUpCheck(disTmp2::Vector{Float64},maxStepTmp::Float64)
disTmp=deepcopy(disTmp2);
dis=disTrans(disTmp);
if findmax(abs,dis)[1]>maxStepTmp
	println("rescale too large displacement")
	disTmp=disTmp./2;
	disUpCheck(disTmp,maxStepTmp);
else
	return disTmp;
end
end

# ╔═╡ fb0006e9-8c3a-4ecd-baaa-8c82024ef2ef
####################in case of too small movement#######################
function disLowCheck(disTmp2::Vector{Float64},minStepTmp::Float64)
disTmp=deepcopy(disTmp2);
dis=disTrans(disTmp);
if findmax(abs,dis)[1]<minStepTmp
	println("rescale too small displacement")
	disTmp=disTmp.*2;
	disLowCheck(disTmp,minStepTmp);
else
	return disTmp;
end
end

# ╔═╡ 354d93d3-b497-40b4-a38d-d9bcc203ffa2
#############in case of too large gap between displacements#############
function reScal(disTmp2::Vector{Float64})
disTmp=deepcopy(disTmp2);
dimension=3;     #for 3D
pNo=Int(length(disTmp)/dimension); #for 3D
dis=disTrans(disTmp);
reDis=Vector{Float64}(undef,pNo);
if pNo==1
	return disTmp
end
minDis=1000000.0;
minInd=1;
for iTmp in 1:pNo
	if minDis>dis[iTmp] && dis[iTmp]!=0
		minDis=dis[iTmp];
		minInd=iTmp;
	end
end
maxDis=findmax(abs,dis)[1];
if maxDis/minDis>50
	for iTmp in 1:pNo
		reDis[iTmp]=sqrt(dis[iTmp]*maxDis);
	end
	for iTmp in 1:pNo
		for jTmp in 1: dimension
			disTmp[(iTmp-1)*dimension+jTmp]=reDis[iTmp]*disTmp[(iTmp-1)*dimension+jTmp]/max(dis[iTmp],0.0000000001);
		end
	end
	reScal(disTmp);
else
	return disTmp
end
end

# ╔═╡ 820e5bfa-b162-4b29-b20b-cb8ad4211348
################in case of ensemble movement#################
function ensMov(disTmp2::Vector{Float64})
disTmp=deepcopy(disTmp2);
countTmp=0;
dimension=3;     #for 3D
pNo=Int(length(disTmp)/dimension);     #for 3D
if pNo==1
	return 0
end
disTmp2=Matrix{Float64}(undef,pNo,dimension);
for iTmp in 1:pNo
	for jTmp in 1:dimension
		disTmp2[iTmp,jTmp]=disTmp[(iTmp-1)*dimension+jTmp];
	end
	if disTmp2[1,:]==disTmp2[iTmp,:]
		countTmp+=1;
	end
end
if countTmp==pNo
	return 1
else
	return 0
end
end

# ╔═╡ 02584a7f-7063-43f2-bf16-144ef5d7562b
###############molecule dynamics for equilibrium################
function md(parPos::Matrix{Float64},Rs::Float64)
parPosTmp=deepcopy(parPos);
precSet=10^-6;	#set the precision of equilibrium position
rndDigNo=9;	#digit No. of round force.
minStepTmp=0.0002;	#minimal step
maxStepTmp=minStepTmp*10;	#maximal step
maxLopNo=500;	#max loop No to find equilibrium position
dimension=3;      #for 3D
ρp=29.0;         #particle density
pNo=length(parPosTmp[:,1]);
parPosData=Array{Float64}(undef,maxLopNo*10,dimension,pNo);
parPosTmp=sepCheck(parPosTmp,Rs);
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
			for iTmp in 1:pNo
				for jTmp in 1:dimension
					parPosTmp[iTmp,jTmp]=round.(parPosTmp[iTmp,jTmp]+disTmp[dimension*(iTmp-1)+jTmp]/2;digits=rndDigNo-2);	#update particles position
					parPosTmp[iTmp,3]=0.0;
					parPosData[countTmp,jTmp,iTmp]=parPosTmp[iTmp,jTmp];
				end
			end
			disTmp2=deepcopy(disTmp);
			break
		else
			for iTmp in 1:pNo
				for jTmp in 1:dimension
					parPosTmp[iTmp,jTmp]=round.(parPosTmp[iTmp,jTmp]+disTmp[dimension*(iTmp-1)+jTmp]/2;digits=rndDigNo-2);	#update particles position
					parPosTmp[iTmp,3]=0.0;
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
return parPosTmp
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
MultipleScattering = "80a8ab25-5750-5d93-a6d7-4adc97cdd5fb"

[compat]
MultipleScattering = "~0.1.17"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.2"
manifest_format = "2.0"
project_hash = "59114702a8d864e98fee96dac8f2160255845187"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "76289dc51920fdc6e0013c872ba9551d54961c24"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.2"
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

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

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

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

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

[[deps.HalfIntegers]]
git-tree-sha1 = "1cfb497b72e1e8ab2256334dee1aaad43fa279ad"
uuid = "f0d1745a-41c9-11e9-1dd9-e5d34d218721"
version = "1.5.1"

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

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.LRUCache]]
git-tree-sha1 = "d36130483e3b6e4cd88d81633b596563264f15db"
uuid = "8ac3fa9e-de4c-5943-b1dc-09c6b5f20637"
version = "1.5.0"

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

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.MultipleScattering]]
deps = ["GSL", "LinearAlgebra", "OffsetArrays", "Printf", "ProgressMeter", "Random", "RecipesBase", "SpecialFunctions", "StaticArrays", "Statistics", "WignerSymbols"]
git-tree-sha1 = "45305f039136d755a786069f08fa04edbf6e0a19"
uuid = "80a8ab25-5750-5d93-a6d7-4adc97cdd5fb"
version = "0.1.17"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "2ac17d29c523ce1cd38e27785a7d23024853a4bb"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.10"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

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
git-tree-sha1 = "00099623ffee15972c16111bcf84c58a0051257c"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.9.0"

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

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

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
deps = ["LinearAlgebra", "Random", "StaticArraysCore"]
git-tree-sha1 = "0adf069a2a490c47273727e029371b31d44b72b2"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.6.5"
weakdeps = ["Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

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

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.WignerSymbols]]
deps = ["HalfIntegers", "LRUCache", "Primes", "RationalRoots"]
git-tree-sha1 = "960e5f708871c1d9a28a7f1dbcaf4e0ee34ee960"
uuid = "9f57e263-0b3d-5e2e-b1be-24f2bb48858b"
version = "2.0.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─6f405d58-bdb6-48ff-8c87-b8eee3a854f5
# ╠═ed3254e0-7b23-11ee-3e45-7dfc3a3015bf
# ╠═05add776-b852-4e6f-9c2f-3ba7726617b0
# ╠═64ddae5f-bea6-4ade-b4b7-ae66f43873f0
# ╠═9c25a47b-e3de-4747-b22a-40a10ded4668
# ╠═fb0006e9-8c3a-4ecd-baaa-8c82024ef2ef
# ╠═354d93d3-b497-40b4-a38d-d9bcc203ffa2
# ╠═820e5bfa-b162-4b29-b20b-cb8ad4211348
# ╠═02584a7f-7063-43f2-bf16-144ef5d7562b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002