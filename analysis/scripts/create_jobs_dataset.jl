using DrWatson

@quickactivate "STLCutters"

using Mustache
using STLCutters
using STLCutters.Tests
s
using STLCutters.Tests: download_things


function jobdict(jobname,params)
  @unpack from, to = params
  includes = 
    "using STLCutters.Tests: run_geometry_list; using FileIO;
    include(\"$(testdir("data","thingi10k_quality_filter.jl"))\");"
  args = "file_ids[$from:$to],load(\"$(tmpdir("thing_to_path.bson"))\")"
  Dict(
  "q" => "normal",
  "o" => datadir(jobname*"o.txt"),
  "e" => datadir(jobname*"e.txt"),
  "walltime" => "24:00:00",
  "ncpus" => 1,
  "mem" => "16gb",
  "name" => jobname,
  "includes" => includes,
  "func" => "run_geometry_list",
  "args" => args,
  "projectdir" => projectdir()
  )
end

template = read(scriptsdir("jobtemplate_gadi.sh"),String)

include(testdir("data","thingi10k_quality_filter.jl")) # file_ids = Int[]

thing_to_path_bson = tmpdir("thing_to_path.bson")

thing_to_path = Dict{Int,Union{Nothing,String}}()
if isfile( thing_to_path_bson )
  thing_to_path = load( thing_to_path_bson )
end

things_ids = file_ids[1:325]
for thing in things_ids
  if !haskey(thing_to_path,thing) || isnothing(thing_to_path[thing])
    filename = download_thing(thing,path=tmpdir())
    thing_to_path[thing] = filename
    save( thing_to_path_bson ,thing_to_path)
  end
end

chunk_size = 100
num_chunks = ceil(Int,length(things_ids)/chunk_size)
starts = [ 1+(i-1)*100 for i in 1:num_chunks ]
ends = [ 100*i for i in 1:num_chunks ]
ends[end] = length(things_ids)

for chunk in 1:num_chunks
  from,to = starts[chunk], ends[chunk]
  range = @dict from to
  jobname = savename(range)
  jobname = replace(jobname,"="=>"_")
  jobfile = datadir(jobname*".sh")
  open(jobfile,"w") do io
    render(io,template,jobdict(jobname,range))
  end
end

