### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ 21b068d0-d091-11ee-142f-65325d95f8fb
"""
run the project:
sweep particle radius, density, soundspeed;
with symmetry.
"""
function sweepParticleMediumProperty()
	#=read basic parameters from input file=#
	input_file=open("input.txt","r");
	output_file=open("output.txt","w")
	readline(input_file);
	readline(input_file);
	δh=parse(Float64,readline(input_file));
	readline(input_file);
	δh2=parse(Float64,readline(input_file));
	readline(input_file);
	δr=parse(Float64,readline(input_file));
	readline(input_file);
	δt=parse(Float64,readline(input_file));
	readline(input_file);
	coef_order1=parse(Int64,readline(input_file));
	readline(input_file);
	coef_order2=parse(Int64,readline(input_file));
	readline(input_file);
	θsteps1=parse(Int64,readline(input_file));
	readline(input_file);
	ϕsteps1=parse(Int64,readline(input_file));
	readline(input_file);
	θsteps2=parse(Int64,readline(input_file));
	readline(input_file);
	ϕsteps2=parse(Int64,readline(input_file));
	readline(input_file);
	max_step_No=parse(Int64,readline(input_file));
	readline(input_file);
	ρb=parse(Float64,readline(input_file));
	readline(input_file);
	cb=parse(Float64,readline(input_file));
	readline(input_file);
	ωb=parse(Float64,readline(input_file));
	ωb=2π*ωb;
	#=read wave parameters=#
	readline(input_file);
	readline(input_file);
	wave_No=parse(Int,readline(input_file));
	inc_amplitude=Array{Float64}(undef,wave_No);
	inc_direction=Matrix{Float64}(undef,wave_No,3);
	inc_position=Matrix{Float64}(undef,wave_No,3);
	for i in 1:wave_No
		readline(input_file);
		inc_amplitude[i]=parse(Float64,readline(input_file));
		readline(input_file);
		direction_str=split(readline(input_file));
		inc_direction[i,:].=[parse(Float64,num) for num in direction_str];
		readline(input_file);
		position_str=split(readline(input_file));
		inc_position[i,:].=[parse(Float64,num) for num in position_str];
		readline(input_file);
	end
	#=read particle parameters=#
	readline(input_file);
	ptc_No=parse(Int64,readline(input_file));
	ptc_position=Matrix{Float64}(undef,ptc_No,3);
	readline(input_file);
	is_center=parse(Bool,readline(input_file));
	readline(input_file);
	ring_No=parse(Int64,readline(input_file));
	readline(input_file);
	radius_str=split(readline(input_file));
	radius_range=[parse(Float64,num) for num in radius_str];
	radius=Array{Float64}(undef,convert(Int64,radius_range[3]));
	for i in 1:convert(Int64,radius_range[3])
		radius[i]=(radius_range[2]-radius_range[1])/radius_range[3]*(i-1)+radius_range[1];
	end
	readline(input_file);
	ρp_str=split(readline(input_file));
	ρp_range=[parse(Float64,num) for num in ρp_str];
	ρp=Array{Float64}(undef,convert(Int64,ρp_range[3]));
	for i in 1:convert(Int64,ρp_range[3])
		ρp[i]=(ρp_range[2]-ρp_range[1])/ρp_range[3]*(i-1)+ρp_range[1];
	end
	readline(input_file);
	cp_str=split(readline(input_file));
	cp_range=[parse(Float64,num) for num in cp_str];
	cp=Array{Float64}(undef,convert(Int64,cp_range[3]));
	for i in 1:convert(Int64,cp_range[3])
		cp[i]=(cp_range[2]-cp_range[1])/cp_range[3]*(i-1)+cp_range[1];
	end
	ring_radius=Array{Float64}(undef,ring_No);
	for i in 1:ring_No
		readline(input_file);
		readline(input_file);
		ring_radius[i]=parse(Float64,readline(input_file));
	end
	write(output_file,"file date information: ",string(today()),"\n");
	writedlm(output_file,["total samples number: ",convert(Int64,radius_range[3]*ρp_range[3]*cp_range[3])]);
	write(output_file,"\n");
	#=calculate force matrix=#
	force_matrix=Matrix{Float64}(undef,2*ptc_No,2*ptc_No);
	generatePosition!(is_center,ptc_No,ring_No,ring_radius,ptc_position);
	for radius_tmp in radius
		radius_tmp2=fill(radius_tmp,ptc_No);
		for ρp_tmp in ρp
			ρp_tmp2=fill(ρp_tmp,ptc_No);
			for cp_tmp in cp
				cp_tmp2=fill(cp_tmp,ptc_No);
				acoustic_model=buildModel(wave_No,ρb,cb,inc_amplitude,inc_direction,inc_position,ptc_No,radius_tmp2,ρp_tmp2,cp_tmp2,ptc_position);
				symmetricStaticMolecularDynamics!(θsteps1,ϕsteps1,ωb,acoustic_model,coef_order1,δr,δh,δt,max_step_No,is_center,ring_No,ring_radius,ptc_position)
				acoustic_model_equilibrium=buildModel2(acoustic_model,ptc_position);
				forceMatrix!(θsteps2,ϕsteps2,ωb,acoustic_model_equilibrium,coef_order2,δr,δh,δh2,force_matrix);
				writedlm(output_file,["particle radius(m):",radius_tmp,"particle density(kg/m^3):",ρp_tmp,"particle soudspeed(m/s):",cp_tmp]);
				write(output_file,"particle equilibrium position(m):\n");
				writedlm(output_file,ptc_position);
				write(output_file,"force matrix(N/m):\n");
				writedlm(output_file,force_matrix,'\t');
				write(output_file,"\n");
				flush(output_file);
			end
		end
	end
	close(input_file);
	close(output_file);
end

# ╔═╡ Cell order:
# ╠═21b068d0-d091-11ee-142f-65325d95f8fb
