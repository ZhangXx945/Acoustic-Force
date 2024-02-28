### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ 90ec2b20-b8d4-11ee-001b-1700b2509f87
"""
generate displacement of each particle from the acoustic force acting on it by static method
"""
function staticParticleDisplacement!(acoustic_model::FrequencySimulation,all_particle_force::Matrix{Float64},all_particle_displacement::Matrix{Float64},δt::Float64)
	ptc_No=length(acoustic_model.particles);
	for i in 1:ptc_No
		for j in 1:2#=for j in 1:3 for 3 dimension displacement=#
			all_particle_displacement[i,j]=3.0/8*all_particle_force[i,j]/π/acoustic_model.particles[i].medium.ρ/(acoustic_model.particles[i].shape.radius)^3*(δt)^2;
		end
		all_particle_displacement[i,3]=0.0;
	end
end

# ╔═╡ 931f665e-f06c-4664-b8bf-b9cd8fcd37cd
"""
generate displacement of each particle from the acoustic force acting on it by dynamic method
"""
function dynamicParticleDisplacement!(acoustic_model::FrequencySimulation,all_particle_force::Matrix{Float64},all_particle_velocity::Matrix{Float64},all_particle_displacement::Matrix{Float64},δt::Float64)
	ptc_No=length(acoustic_model.particles);
	for i in 1:ptc_No
		for j in 1:2#=for j in 1:3 for 3 dimension displacement=#
			all_particle_displacement[i,j]=all_particle_velocity[i,j]*δt+3.0/8*all_particle_force[i,j]/π/acoustic_model.particles[i].medium.ρ/(acoustic_model.particles[i].shape.radius)^3*(δt)^2;
			all_particle_velocity[i,j]=all_particle_velocity[i,j]+3.0/4*all_particle_force[i,j]/π/acoustic_model.particles[i].medium.ρ/(acoustic_model.particles[i].shape.radius)^3*δt;
		end
		all_particle_displacement[i,3]=0.0;
		all_particle_velocity[i,3]=0.0;
	end
end

# ╔═╡ 0bea4125-c36b-49b2-8ed0-0a4edccd2bd4
"""
check if any two particles collide in the system, and reorganize their positions if any overlapping occurs
"""
function checkOverlap!(acoustic_model::FrequencySimulation,ptc_position::Matrix{Float64},δr)
	ptc_No=length(acoustic_model.particles);
	if ptc_No==1
		return Nothing
	end
	#=the following code is for 3 dimension case
	for i in 1:ptc_No-1
		for j in (i+1):ptc_No
			sep_distance=sqrt((ptc_position[j,1]-ptc_position[i,1])^2+(ptc_position[j,2]-ptc_position[i,2])^2+(ptc_position[j,3]-ptc_position[i,3])^2);
			sep_minimum=acoustic_model.particles[i].shape.radius+acoustic_model.particles[j].shape.radius+δr
			if sep_dis<sep_minimum
				θ=atan(sqrt((ptc_position[j,1]-ptc_position[i,1])^2+(ptc_position[j,2]-ptc_position[i,2])^2),ptc_position[j,3]-ptc_position[i,3]);
				ϕ=atan(ptc_position[j,2]-ptc_position[i,2],ptc_position[j,1]-ptc_position[i,1]);
				ptc_position[j,1]=2*sep_minimum*sin(θ)*cos(ϕ);
				ptc_position[j,2]=2*sep_minimum*sin(θ)*sin(ϕ);
				ptc_position[j,3]=2*sep_minimum*cos(θ);=#
	for i in 1:ptc_No-1
		for j in (i+1):ptc_No
			sep_distance=sqrt((ptc_position[j,1]-ptc_position[i,1])^2+(ptc_position[j,2]-ptc_position[i,2])^2);
			sep_minimum=acoustic_model.particles[i].shape.radius+acoustic_model.particles[j].shape.radius+δr
			if sep_distance<sep_minimum
				println("Overlapping between particle ",i," and ",j)
				ϕ=atan(ptc_position[j,2]-ptc_position[i,2],ptc_position[j,1]-ptc_position[i,1]);
				ptc_position[j,1]=2*sep_minimum*cos(ϕ);
				ptc_position[j,2]=2*sep_minimum*sin(ϕ);
			end
		end
	end
end 

# ╔═╡ 74917bc8-8d84-45d4-8f41-cbd65b7f2e9b
"""
find the minimal non-zero number from an array
"""
function minNonzeroNumber(arr::Array{Float64})
	nonzero_indices = findall(x -> x != 0, arr)
	if isempty(nonzero_indices)
        return 0.0
    end
	min_nonzero=findmin(arr[nonzero_indices])[1];
	return min_nonzero
end

# ╔═╡ 94432750-86e1-4604-b1a1-b514fa7f02a9
"""
rescale displacement if necessary
"""
function rescaleDisplacement!(all_particle_displacement::Matrix{Float64},max_displacement::Float64,min_displacement::Float64,displacement_gap::Float64)
	max_displacement_tmp=findmax(sqrt.(sum((all_particle_displacement.*all_particle_displacement),dims=2)))[1];
	min_nonzero_displacement_tmp=minNonzeroNumber(sqrt.(sum((all_particle_displacement.*all_particle_displacement),dims=2)));
	ptc_No=size(all_particle_displacement)[1];
	#=in case of too large displacement=#
	if max_displacement_tmp>max_displacement
		all_particle_displacement./=2;
		rescaleDisplacement!(all_particle_displacement,max_displacement,min_displacement,displacement_gap);
	#=in case of too small displacement=#
	elseif max_displacement_tmp<min_displacement
		all_particle_displacement.*=2;
		rescaleDisplacement!(all_particle_displacement,max_displacement,min_displacement,displacement_gap);
	end
	if ptc_No==1
		return Nothing
	end
	if min_nonzero_displacement_tmp==0.0
		return Nothing;
	end
	if 	max_displacement_tmp/min_nonzero_displacement_tmp>displacement_gap
		for i in 1:ptc_No
			#=for two dimensional movement only=#
			ϕ=atan(all_particle_displacement[i,2],all_particle_displacement[i,1]);
			all_particle_displacement[i,1]=cbrt.(sqrt(sum(all_particle_displacement[i,:].*all_particle_displacement[i,:])))*cos(ϕ);
			all_particle_displacement[i,2]=cbrt.(sqrt(sum(all_particle_displacement[i,:].*all_particle_displacement[i,:])))*sin(ϕ);
			all_particle_displacement[i,3]=0.0;
		end
		max_displacement_tmp=findmax(sqrt.(sum((all_particle_displacement.*all_particle_displacement),dims=2)))[1];
		min_nonzero_displacement_tmp=minNonzeroNumber(sqrt.(sum((all_particle_displacement.*all_particle_displacement),dims=2)));
		rescaleDisplacement!(all_particle_displacement,max_displacement,min_displacement,displacement_gap);
	else
		return Nothing
	end
end		

# ╔═╡ b6845361-554d-4328-b048-94a22c86de7b
"""
in case of ensemble motion of particles
"""
function ensembleMotion(all_particle_displacement::Matrix{Float64})
	ptc_No=size(all_particle_displacement)[1];
	count=1;
	for i in 1:ptc_No-1
		if all_particle_displacement[i,1]==all_particle_displacement[i+1,1] && all_particle_displacement[i,2]==all_particle_displacement[i+1,2]
			count=count+1;
		end
	end
	if count==ptc_No
		println("EQUILIBRIUM, ENSEMBLE MOVEMENT")
		return 1
	end
end

# ╔═╡ 9ba7f962-a70c-48ac-ab8f-2943c773894e
"""
static molecular dynamics driven by acoustic radiation force
"""
function staticMolecularDynamics!(θsteps::Int64,ϕsteps::Int64,ωb::Float64,acoustic_model::FrequencySimulation,coef_order::Int64,δr::Float64,δh::Float64,δt::Float64,max_step_No::Int64,ptc_position::Matrix{Float64})
	#======================initialization===================#
	ptc_No=length(acoustic_model.particles);
	ptc_force=Matrix{Float64}(undef,ptc_No,3);#=store force data=#
	ptc_disp1=Matrix{Float64}(undef,ptc_No,3);#=store displacement data=#
	ptc_disp2=deepcopy(ptc_disp1);
	radius=Array{Float64}(undef,ptc_No);#=store particles radius=#
	λb=2π*real(acoustic_model.source.medium.c)/ωb;
	ptc_traj1=Array{Float64}(undef,ptc_No,3,0);#=store particles trajectory=#
	ptc_traj2=deepcopy(ptc_traj1);
	for i in 1:ptc_No
		radius[i]=acoustic_model.particles[i].shape.radius;
	end
	up_stepsize=findmin([findmin(radius)[1],λb])[1]/40;
	low_stepsize=up_stepsize/10;
	min_stepsize=low_stepsize/100;
	step_count=0;
	round_count=0;
	trunc_digits=1;
	while convert(Float64,1/(10^trunc_digits))>min_stepsize
		trunc_digits+=1;
	end
	trunc_digits+=5;
	checkOverlap!(acoustic_model,ptc_position,δr);
	ptc_position.=trunc.(ptc_position,digits=trunc_digits);
	#==================molecular dynamics=========================#
	while low_stepsize>=min_stepsize
		#=========each round==========#
		for i in 1:max_step_No
			#===========each step========#
			acoustic_model=buildModel2(acoustic_model,ptc_position);
			scattering_coefficient=getCoefficients(acoustic_model,ωb,coef_order);
			allParticleAcousticForce!(θsteps,ϕsteps,ωb,acoustic_model,scattering_coefficient,δr,δh,ptc_force);
			if step_count!=0
				ptc_disp2.=ptc_disp1;
			end 
			staticParticleDisplacement!(acoustic_model,ptc_force,ptc_disp1,δt);
			ptc_disp1.=trunc.(ptc_disp1,digits=trunc_digits);
			if findmax(abs,ptc_disp1)[1]==0.0
				println("particles are in EQUILIBRIUM positions!")
				return Nothing
			end
			if ensembleMotion(ptc_disp1)==1
				println("ensemble motion begins!")
				return Nothing
			end
			rescaleDisplacement!(ptc_disp1,up_stepsize,low_stepsize,100.0);
			ptc_position.+=ptc_disp1;
			ptc_position.=trunc.(ptc_position,digits=trunc_digits);
			checkOverlap!(acoustic_model,ptc_position,δr);
			if step_count!=0
				ptc_traj2=deepcopy(ptc_traj1);
			end
			step_count+=1;
			ptc_traj1=Array{Float64}(undef,ptc_No,3,step_count);
			for i in 1:step_count-1
				ptc_traj1[:,:,i].=ptc_traj2[:,:,i];
			end
			ptc_traj1[:,:,step_count].=ptc_position;
			#===================output information=======================#
			println("Max Stepsize: ",up_stepsize,"; Total Step No: ",step_count);
			println("Round No: ",round_count,"; Current Round Step No: ",i);
			println("Displacement:\n",ptc_disp1,";\n Position:\n",ptc_position);
			#====================plot 3D diagram==========================#
			plot3d(plottitles=step_count);
			for j in 1:ptc_No
			display(path3d!(ptc_traj1[j,1,:],ptc_traj1[j,2,:],ptc_traj1[j,3,:]))
			end
			if i==max_step_No
				println("reach MAXIMAL step number!")
			end
			#===================break current loop======================#
			if step_count!=0 && findmax(ptc_disp1.*ptc_disp2)[1]<=0
				println("particles are in EQUILIBRIUM positions!")
				break
			end
		end
		round_count+=1;
		up_stepsize=up_stepsize/2;
		low_stepsize=low_stepsize/2;
	end		
end	

# ╔═╡ db81b7be-0907-4745-86e8-932903b3ddd2
"""
generate particle positions from symmetric information
"""
function generatePosition!(is_center::Bool,ptc_No::Int64,ring_No::Int64,ring_radius::Array{Float64},ptc_position::Matrix{Float64})
	if is_center
		ptc_position[1,1]=0.0;
		ptc_position[1,2]=0.0;
		ptc_position[1,3]=0.0;
		if (ptc_No-1)%ring_No==0
			ring_ptc_No=convert(Int64,(ptc_No-1)/ring_No);
		else
			println("particles cannot be evenly distributed in rings!")
			return Nothing	
		end
		for i in 1:ring_No
			for j in 1:ring_ptc_No
				ptc_position[1+j+(i-1)*ring_ptc_No,1]=ring_radius[i]*cos(2π/ring_ptc_No*(j-1));
				ptc_position[1+j+(i-1)*ring_ptc_No,2]=ring_radius[i]*sin(2π/ring_ptc_No*(j-1));
				ptc_position[1+j+(i-1)*ring_ptc_No,3]=0.0;
			end
		end
	else
		if ptc_No%ring_No==0
			ring_ptc_No=convert(Int64,ptc_No/ring_No);
		else
			println("particles cannot be evenly distributed in rings!")
			return Nothing
		end
		for i in 1:ring_No
			for j in 1:ring_ptc_No
				ptc_position[j+(i-1)*ring_ptc_No,1]=ring_radius[i]*cos(2π/ring_ptc_No*(j-1));
				ptc_position[j+(i-1)*ring_ptc_No,2]=ring_radius[i]*sin(2π/ring_ptc_No*(j-1));
				ptc_position[j+(i-1)*ring_ptc_No,3]=0.0;
			end
		end
	end
end

# ╔═╡ d7513a9c-9f7e-4031-a5cb-f1c5968daf4d
"""
static molecular dynamics driven by acoustic radiation force for symmetric system
"""
function symmetricStaticMolecularDynamics!(θsteps::Int64,ϕsteps::Int64,ωb::Float64,acoustic_model::FrequencySimulation,coef_order::Int64,δr::Float64,δh::Float64,δt::Float64,max_step_No::Int64,is_center::Bool,ring_No::Int64,ring_radius::Array{Float64},ptc_position::Matrix{Float64})

	#===================initialization================#
	ptc_No=length(acoustic_model.particles);
	generatePosition!(is_center,ptc_No,ring_No,ring_radius,ptc_position);
	if is_center
		if (ptc_No-1)%ring_No==0
			ring_ptc_No=convert(Int64,(ptc_No-1)/ring_No);
		else
			println("particles cannot be evenly distributed in rings!")
			return Nothing
		end
	else
		if ptc_No%ring_No==0
			ring_ptc_No=convert(Int64,ptc_No/ring_No);
		else
			println("particles cannot be evenly distributed in rings!")
			return Nothing
		end
	end
	ptc_force=Matrix{Float64}(undef,ptc_No,3);#=store force data=#
	ptc_disp1=Matrix{Float64}(undef,ptc_No,3);#=store displacement data=#
	ptc_disp2=deepcopy(ptc_disp1);
	radius=Array{Float64}(undef,ptc_No);#=store particles radius=#
	λb=2π*real(acoustic_model.source.medium.c)/ωb;
	ptc_traj1=Array{Float64}(undef,ptc_No,3,0);#=store particles trajectory=#
	ptc_traj2=deepcopy(ptc_traj1);
	for i in 1:ptc_No
		radius[i]=acoustic_model.particles[i].shape.radius;
	end
	up_stepsize=findmin([findmin(radius)[1],λb])[1]/40;
	low_stepsize=up_stepsize/10;
	min_stepsize=low_stepsize/100;
	step_count=0;
	round_count=1;
	trunc_digits=1;
	while convert(Float64,1/(10^trunc_digits))>min_stepsize
		trunc_digits+=1;
	end
	trunc_digits+=5;
	#=================in case of overlapping=================#
	if is_center
		if ring_radius[1]<=(radius[2]+radius[1]+2δr)
			println("OVERLAPPING")
			ring_radius[1]=2*(radius[2]+radius[1]+2δr)
		end
		for i in 1:(ring_No-1)
			if ptc_position[2+i*ring_ptc_No,1]-ptc_position[2+(i-1)*ring_ptc_No,1]<=radius[2+i*ring_ptc_No]+radius[2+(i-1)*ring_ptc_No]+2*δr
				println("OVERLAPPING")
				ring_radius[i+1]=2*(radius[2+i*ring_ptc_No]+radius[2+(i-1)*ring_ptc_No]+2*δr);
				generatePosition!(is_center,ptc_No,ring_No,ring_radius,ptc_position);
			end
		end
	else
		if ptc_position[1,1]<=((radius[1]+δr)/sin(π/ring_ptc_No))
			println("OVERLAPPING")
			ring_radius[1]=2*(radius[1]+δr)/sin(π/ring_ptc_No);
			generatePosition!(is_center,ptc_No,ring_No,ring_radius,ptc_position);
		end
		for i in 1:(ring_No-1)
			if ptc_position[1+i*ring_ptc_No,1]-ptc_position[1+(i-1)*ring_ptc_No,1]<=radius[1+i*ring_ptc_No]+radius[1+(i-1)*ring_ptc_No]+2*δr
				println("OVERLAPPING")
				ring_radius[i+1]=2*(radius[1+i*ring_ptc_No]+radius[1+(i-1)*ring_ptc_No]+2*δr);
				generatePosition!(is_center,ptc_No,ring_No,ring_radius,ptc_position);
			end
		end
	end
	ptc_position.=trunc.(ptc_position,digits=trunc_digits);
	#==================molecular dynamics=========================#
	while low_stepsize>=min_stepsize
		#=========each round==========#
		for i in 1:max_step_No
			#===========each step========#
			acoustic_model=buildModel2(acoustic_model,ptc_position);
			scattering_coefficient=getCoefficients(acoustic_model,ωb,coef_order);
			allParticleAcousticForce!(θsteps,ϕsteps,ωb,acoustic_model,scattering_coefficient,δr,δh,ptc_force);
			if step_count!=0
				ptc_disp2.=ptc_disp1;
			end 
			staticParticleDisplacement!(acoustic_model,ptc_force,ptc_disp1,δt);
			ptc_disp1.=trunc.(ptc_disp1,digits=trunc_digits);
			if findmax(abs,ptc_disp1)[1]==0.0
				println("particles are in EQUILIBRIUM positions!")
				return Nothing
			end
			if ensembleMotion(ptc_disp1)==1
				println("ensemble motion begins!")
				return Nothing
			end
			rescaleDisplacement!(ptc_disp1,up_stepsize,low_stepsize,100.0);
			ptc_position.+=ptc_disp1;
			if is_center
				for i in 1:ring_No
					ring_radius[i]=ptc_position[2+ring_ptc_No*(i-1),1]
				end
			else
				for i in 1:ring_No
					ring_radius[i]=ptc_position[1+ring_ptc_No*(i-1),1]
				end
			end
			generatePosition!(is_center,ptc_No,ring_No,ring_radius,ptc_position);
			ptc_position.=trunc.(ptc_position,digits=trunc_digits);
			#=================in case of overlapping=================#
	if is_center
		if ring_radius[1]<=(radius[2]+radius[1]+2δr)
			println("OVERLAPPING")
			ring_radius[1]=2*(radius[2]+radius[1]+2δr)
		end
		for i in 1:(ring_No-1)
			if ptc_position[2+i*ring_ptc_No,1]-ptc_position[2+(i-1)*ring_ptc_No,1]<=radius[2+i*ring_ptc_No]+radius[2+(i-1)*ring_ptc_No]+2*δr
				println("OVERLAPPING")
				ring_radius[i+1]=2*(radius[2+i*ring_ptc_No]+radius[2+(i-1)*ring_ptc_No]+2*δr);
				generatePosition!(is_center,ptc_No,ring_No,ring_radius,ptc_position);
			end
		end
	else
		if ptc_position[1,1]<=((radius[1]+δr)/sin(π/ring_ptc_No))
			println("OVERLAPPING")
			ring_radius[1]=2*(radius[1]+δr)/sin(π/ring_ptc_No);
			generatePosition!(is_center,ptc_No,ring_No,ring_radius,ptc_position);
		end
		for i in 1:(ring_No-1)
			if ptc_position[1+i*ring_ptc_No,1]-ptc_position[1+(i-1)*ring_ptc_No,1]<=radius[1+i*ring_ptc_No]+radius[1+(i-1)*ring_ptc_No]+2*δr
				println("OVERLAPPING")
				ring_radius[i+1]=2*(radius[1+i*ring_ptc_No]+radius[1+(i-1)*ring_ptc_No]+2*δr);
				generatePosition!(is_center,ptc_No,ring_No,ring_radius,ptc_position);
			end
		end
	end
	ptc_position.=trunc.(ptc_position,digits=trunc_digits);
			#=================in case of overlapping=================#
			if step_count!=0
				ptc_traj2=deepcopy(ptc_traj1);
			end
			step_count+=1;
			ptc_traj1=Array{Float64}(undef,ptc_No,3,step_count);
			for i in 1:step_count-1
				ptc_traj1[:,:,i].=ptc_traj2[:,:,i];
			end
			ptc_traj1[:,:,step_count].=ptc_position;
			#===================output information=======================#
			println("Max Stepsize: ",up_stepsize,"; Total Step No: ",step_count);
			println("Round No: ",round_count,"; Current Round Step No: ",i);
			println("Displacement:\n",ptc_disp1,";\nPosition:\n",ptc_position);
			#====================plot 3D diagram==========================#
			plot3d(plottitles=step_count);
			for j in 1:ptc_No
			display(path3d!(ptc_traj1[j,1,:],ptc_traj1[j,2,:],ptc_traj1[j,3,:]))
			end
			if i==max_step_No
				println("reach MAXIMAL step number!")
			end
			#===================break current loop======================#
			if step_count!=0 && findmax(ptc_disp1.*ptc_disp2)[1]<=0.0
				println("particles are in EQUILIBRIUM positions!")
				break
			end
		end
		round_count+=1;
		up_stepsize=up_stepsize/2;
		low_stepsize=low_stepsize/2;
	end		
end

# ╔═╡ Cell order:
# ╠═90ec2b20-b8d4-11ee-001b-1700b2509f87
# ╠═931f665e-f06c-4664-b8bf-b9cd8fcd37cd
# ╠═0bea4125-c36b-49b2-8ed0-0a4edccd2bd4
# ╠═74917bc8-8d84-45d4-8f41-cbd65b7f2e9b
# ╠═94432750-86e1-4604-b1a1-b514fa7f02a9
# ╠═b6845361-554d-4328-b048-94a22c86de7b
# ╠═9ba7f962-a70c-48ac-ab8f-2943c773894e
# ╠═db81b7be-0907-4745-86e8-932903b3ddd2
# ╠═d7513a9c-9f7e-4031-a5cb-f1c5968daf4d
