module CreateJobsDataset

using DrWatson

@quickactivate "STLCuttersAnalysis"

using Mustache
using STLCuttersAnalysis

thingi10_solid_manifold = datadir("geometries","thingi10k_solid_manifold.jl")

include( thingi10_solid_manifold ) # file_ids::Vector{Int}

include( scriptsdir("generic","hpc_dicts.jl") )

function testdict(hpcname,from,to)
  incl = 
      "using STLCutters.Tests: run_geometry_list; using FileIO;
       include(\"$thingi10_solid_manifold\");"
  func = "run_geometry_list"
  args = "
        file_ids[$from:$to],
        load(\"$(tmpdir("thing_to_path.bson"))\"),"
  kwargs = "
        rerun = false,
        datapath = \"$(datadir(hpcname))\""
  @dict incl func args kwargs
end

function jobdir(params)
  jobname = savename(params)
  jobname = replace(jobname,"="=>"_")
end

function main(
  hpc_id,
  subset;
  walltime="24:00:00",
  memory=16,
  ncpus=memory÷4,
  chunk_size=100)

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

  num_chunks = ceil(Int,length(things_ids)/chunk_size)
  starts = [ 1+(i-1)*chunk_size for i in 1:num_chunks ]
  ends = [ chunk_size*i for i in 1:num_chunks ]
  ends[end] = length(things_ids)

  jobparams = @dict walltime memory ncpus
  for chunk in 1:num_chunks
    from,to = starts[chunk], ends[chunk]
    from,to = subset[from], subset[to]
    range = @dict from to
    jobname = jobdir(range)
    jobfile = datadir(hpcname,jobname*".sh")
    @pack! jobparams = jobname
    testparams = testdict(hpcname,from,to)
    open(jobfile,"w") do io
      render(io,template,jobdict(hpc_id,jobparams,testparams))
    end
  end
end

end # module

