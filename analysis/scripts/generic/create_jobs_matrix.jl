module CreateJobsMatrix

using DrWatson

@quickactivate "STLCutters"

using Mustache
using STLCutters
using STLCutters.Tests

include( scriptsdir("generic","hpc_dicts.jl") )

function jobdict(hpc_id,jobname,params)
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
  "walltime" => "24:00:00",
  "ncpus" => 1,
  "mem" => memory(16),
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
  "441708", # Standford Bunny
  "35269", # Octocat
  "65904", # Heart 
  "551021", # Art de Triomph
  "37266", # Extruded Earth
  "252119"] # Angel

function create_jobs_matrix(hpc_id)

  @unpack hpcname = hpc_dict[hpc_id]

  template = read(scriptsdir(hpcname,"jobtemplate.sh"),String)

  kwargs = 
         "displacements=exp10.(-17:-1 ),
          angles=exp10.(-17:-1 ),
          rerun=false,
          datapath = \"$(datadir(hpcname))\""

  all_params = Dict(
    :func => "rotations_and_displacements",
    :filename => stl_list,
    :nmin => 1,
    :nmax => 14 .* 2 .^ (0:5),
    :kwargs => kwargs 
    )

  if !isdir( datadir(hpcname) )
    mkdir( datadir(hpcname) )
  end

  for params in dict_list(all_params)
    jobname = savename(params[:filename],params,ignores=(:func,:filename,:kwargs))
    jobname = replace(jobname,"="=>"")
    jobfile = datadir(hpcname,jobname*".sh")
    open(jobfile,"w") do io
      render(io,template,jobdict(hpc_id,jobname,params))
    end
  end

end

export create_jobs_matrix

end # module

using Main.CreateJobsMatrix 
