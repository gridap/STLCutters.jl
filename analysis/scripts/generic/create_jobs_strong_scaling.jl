using DrWatson

@quickactivate "STLCutters"

using Mustache
using STLCutters
using STLCutters.Tests

include( scriptsdir("generic","hpc_dicts.jl") )

function jobdict(hpc_id,jobname,params,walltime,mem_gb)
  @unpack func,filename,nmin,nmax,ncpus,threading,kwargs = params
  @unpack hpcname,queue,memory,julia= hpc_dict[hpc_id]
  filename = testdir("data","$filename.stl")
  args = "
          \"$filename\",
          nmin=$nmin,
          nmax=$nmax,
          threading=:$threading,
          nthreads=$ncpus,
          $kwargs"
  mem_gb = max(mem_gb,ncpus*4)
  Dict(
  "q" => queue,
  "o" => datadir(hpcname,jobname*".out"),
  "e" => datadir(hpcname,jobname*".err"),
  "walltime" => walltime,
  "ncpus" => ncpus,
  "mem" => "190gb",#memory(mem_gb),
  "julia" => julia,
  "name" => jobname,
  "func" => func,
  "args" => args,
  "projectdir" => projectdir()
  )
end

function create_jobs_strong_scaling(
  hpc_id;
  walltime="1:00:00",
  memory=16,
  ncpus=[1,2,4,8,16,32,48],
  nmaxs=100,
  stls)

  @unpack hpcname = hpc_dict[hpc_id]

  template = read(scriptsdir(hpcname,"jobtemplate.sh"),String)

  
  kwargs = 
         "nruns=5,
          rerun=false,
          datapath = \"$(datadir(hpcname))\""

  func = "run_and_save"

  prefix = "strong_"

  all_params = Dict(
    :func => func,
    :filename => stls,
    :nmin => 1,
    :nmax => nmaxs,
    :ncpus => ncpus,
    :threading => [:threads,:spawn],
    :kwargs => kwargs 
    )
 
  for params in dict_list(all_params)
    jobname = savename(prefix*params[:filename],params,ignores=(:func,:filename,:kwargs))
    jobname = replace(jobname,"="=>"")
    jobname = replace(jobname,"threading"=>"")
    jobfile = datadir(hpcname,jobname*".sh")
    open(jobfile,"w") do io
      render(io,template,jobdict(hpc_id,jobname,params,walltime,memory))
    end
  end
end
   
## cube,550964,35269,441708
## create_jobs_strong_scaling(:gadi,stls=["cube,441708"],nmaxs=[100,200])
