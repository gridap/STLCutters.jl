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
  incl = ""
  kwargs = "
          nmin=$nmin,
          nmax=$nmax,
          $kwargs"
  args = "\"$filename\",\n"
  testparams = Dict{Symbol,Any}()
  @dict incl func args kwargs
end

function jobdir(prefix,params)
  ignores = (:func,:filename,:kwargs)
  jobname = savename(prefix*params[:filename],params;ignores)
  replace(jobname,"="=>"")
end

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

  _kwargs = join(["$k = $v, " for (k,v) in kwargs ])

  _kwargs = 
         "$_kwargs
          rerun=false,
          datapath = \"$(datadir(hpcname))\","

  if displace
    func = "rotations_and_displacements"
    _kwargs =
         "$_kwargs
          displacements=exp10.(-17:-1 ),
          angles=exp10.(-17:-1 ),"
  else
    func = "run_and_save"
  end

  if poisson
    _kwargs =
         "$_kwargs
          poisson=true,
          agfem_threshold=0.5,"
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
  push!(all_params,kwargs...)

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

