module CreateJobsMatrix

using DrWatson

@quickactivate "STLCuttersAnalysis"

using Mustache

include( projectdir("scripts","generic","hpc_dicts.jl") )

stl_list = [
  "cube",
  "37881", # Gear
  "47076", # Azteca pyramid
  "96457", # NightHawk
  "550964", # Eiffel Tower
  "293137", # Bunny low-poly
  "441708", # Standford Bunny
  "35269", # Octocat
  "65904_1", # Heart 
  "551021", # Arc de Triomph
  "37266", # Extruded Earth
  "252119"] # Angel

function testdict(params)
  @unpack func,filename,nmin,nmax,kwargs = params 
  filename = datadir("geometries","$filename.stl")
  incl = ""
  kwargs = 
   s10 * "nmin=$nmin,\n" *
   s10 * "nmax=$nmax,\n" *
   kwargs
  args = "\"$filename\",\n"
  testparams = Dict{Symbol,Any}()
  @dict incl func args kwargs
end

function jobdir(prefix,params)
  ignores = (:func,:filename,:kwargs)
  jobname = savename(prefix*params[:filename],params;ignores)
  replace(jobname,"="=>"")
end

format_data(a::String) = "\"$a\""

format_data(a::Number) = "$a"

format_data(a::Symbol) = ":$a"

s10 = join(fill(' ',10))

function main(
  hpc_id;
  walltime="24:00:00",
  memory=16,
  ncpus=4,
  nmaxs=14 .* 2 .^ (0:5),
  displace=true,
  poisson=false,
  kwargs...)

  @unpack hpcname = hpc_dict[hpc_id]

  template = read( projectdir("scripts",hpcname,"jobtemplate.sh"),String)

  _kwargs = join([ s10*"$k=$(format_data(v)),\n" for (k,v) in kwargs ])

  _kwargs *= s10*"rerun=false,\n"
  _kwargs *= s10*"datapath = \"$(datadir(hpcname))\",\n"

  if displace
    func = "rotations_and_displacements"
    _kwargs *= s10*"displacements=exp10.(-17:-1 ),\n"
    _kwargs *= s10*"angles=exp10.(-17:-1 ),\n"
  else
    func = "run_and_save"
  end

  if poisson
    _kwargs *= s10*"poisson=true,\n"
    _kwargs *= s10*"agfem_threshold=0.5,\n"
  end

  prefix = ""
  if poisson
    prefix = "poisson_"
  end
  if !displace
    prefix = prefix*"single_"
  end

  all_params = Dict(
    :func => func,
    :filename => stl_list,
    :nmin => 1,
    :nmax => nmaxs,
    :kwargs => _kwargs 
    )
  for kwarg in kwargs
    push!(all_params,kwarg)
  end

  if !isdir( datadir(hpcname) )
    mkdir( datadir(hpcname) )
  end

  jobparams = @dict walltime memory ncpus
  for params in dict_list(all_params)
    jobname = jobdir(prefix,params)
    jobfile = datadir(hpcname,jobname*".sh")
    @pack! jobparams = jobname
    testparams = testdict(params)
    open(jobfile,"w") do io
      render(io,template,jobdict(hpc_id,jobparams,testparams))
    end
  end
end

end # module

