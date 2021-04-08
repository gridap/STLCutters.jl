module CreateJobsMatrix

using DrWatson

@quickactivate "STLCutters"

using Mustache
using STLCutters
using STLCutters.Tests

include( scriptsdir("generic","hpc_dicts.jl") )

function jobdict(hpc_id,jobname,params,walltime,mem_gb,ncpus)
  @unpack func,filename,nmin,nmax,kwargs = params
  @unpack hpcname,queue,memory,julia= hpc_dict[hpc_id]
  filename = testdir("data","$filename.stl")
  args = "
          \"$filename\",
          nmin=$nmin,
          nmax=$nmax,
          $kwargs"
  Dict(
  "q" => queue,
  "o" => datadir(hpcname,jobname*".out"),
  "e" => datadir(hpcname,jobname*".err"),
  "walltime" => walltime,
  "ncpus" => ncpus,
  "mem" => memory(mem_gb),
  "julia" => julia,
  "name" => jobname,
  "func" => func,
  "args" => args,
  "projectdir" => projectdir()
  )
end

stl_list = [
  "cube",
  "wine_glass",
  "37881", # Gear
  "47076", # Azteca pyramid
  "96457", # NightHawk
  "550964", # Eiffel Tower
  "293137", # Bunny low-poly
  "441708_1", # Standford Bunny
  "35269", # Octocat
  "65904", # Heart 
  "551021", # Arc de Triomph
  "37266", # Extruded Earth
  "252119"] # Angel

function create_jobs_matrix(
  hpc_id;
  walltime="24:00:00",
  memory=16,
  ncpus=4,
  nmaxs=14 .* 2 .^ (0:5),
  displace=true,
  poisson=false,
  solution_order=1,
  solver=:direct)

  @unpack hpcname = hpc_dict[hpc_id]

  template = read(scriptsdir(hpcname,"jobtemplate.sh"),String)

  kwargs = 
         "rerun=false,
          datapath = \"$(datadir(hpcname))\""

  if displace
    func = "rotations_and_displacements"
    kwargs =
         "$kwargs,
          displacements=exp10.(-17:-1 ),
          angles=exp10.(-17:-1 )"
  else
    func = "run_and_save"
  end

  if poisson
    kwargs =
         "$kwargs,
          poisson=true,
          solution_order=$solution_order,
          agfem_threshold=0.5,
          solver=$solver"
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
    :kwargs => kwargs 
    )

  if !isdir( datadir(hpcname) )
    mkdir( datadir(hpcname) )
  end

  for params in dict_list(all_params)
    jobname = savename(prefix*params[:filename],params,ignores=(:func,:filename,:kwargs))
    jobname = replace(jobname,"="=>"")
    jobfile = datadir(hpcname,jobname*".sh")
    open(jobfile,"w") do io
      render(io,template,jobdict(hpc_id,jobname,params,walltime,memory,ncpus))
    end
  end

end

export create_jobs_matrix

end # module

import Main.CreateJobsMatrix: create_jobs_matrix
