module CreateJobsDataset

using DrWatson

@quickactivate "STLCutters"

using Mustache
using STLCutters
using STLCutters.Tests

using STLCutters.Tests: download_things

export create_jobs_dataset

function jobdict(hpc_id,jobname,params,walltime,mem_gb)
  @unpack from, to = params
  @unpack hpcname,queue,memory,julia= hpc_dict[hpc_id]
  includes = 
      "using STLCutters.Tests: run_geometry_list; using FileIO;
       include(\"$(testdir("data","thingi10k_quality_filter.jl"))\");"
  func = "run_geometry_list"
  args = "
        file_ids[$from:$to],
        load(\"$(tmpdir("thing_to_path.bson"))\"),
        rerun = false,
        datapath = \"$(datadir(hpcname))\""
  Dict(
  "q" => queue,
  "o" => datadir(hpcname,jobname*".out"),
  "e" => datadir(hpcname,jobname*".err"),
  "walltime" => walltime,
  "ncpus" => 1,
  "mem" => memory(mem_gb),
  "name" => jobname,
  "julia" => julia,
  "includes" => includes,
  "func" => func,
  "args" => args,
  "projectdir" => projectdir()
  )
end

include( testdir("data","thingi10k_quality_filter.jl") ) # file_ids = Int[]

include( scriptsdir("generic","hpc_dicts.jl") )

function create_jobs_dataset(hpc_id,subset;walltime="24:00:00",memory=16)

  @unpack hpcname = hpc_dict[hpc_id]

  template = read(scriptsdir(hpcname,"jobtemplate.sh"),String)


  thing_to_path_bson = tmpdir("thing_to_path.bson")

  thing_to_path = Dict{Int,Union{Nothing,String}}()
  if isfile( thing_to_path_bson )
    thing_to_path = load( thing_to_path_bson )
  end

  things_ids = file_ids[ subset ]
  for thing in things_ids
    if !haskey(thing_to_path,thing) || isnothing(thing_to_path[thing])
      filename = download_thing(thing,path=tmpdir())
      thing_to_path[thing] = filename
      save( thing_to_path_bson ,thing_to_path)
    end
  end

  chunk_size = 100
  num_chunks = ceil(Int,length(things_ids)/chunk_size)
  starts = [ 1+(i-1)*chunk_size for i in 1:num_chunks ]
  ends = [ chunk_size*i for i in 1:num_chunks ]
  ends[end] = length(things_ids)

  for chunk in 1:num_chunks
    from,to = starts[chunk], ends[chunk]
    from,to = subset[from], subset[to]
    range = @dict from to
    jobname = savename(range)
    jobname = replace(jobname,"="=>"_")
    jobfile = datadir(hpcname,jobname*".sh")
    open(jobfile,"w") do io
      render(io,template,jobdict(hpc_id,jobname,range,walltime,memory))
    end
  end

end

end # module

using Main.CreateJobsDataset

