using DrWatson

@quickactivate "STLCutters"

using Mustache
using STLCutters
using STLCutters.Tests

function jobdict(jobname,params)
  filename = testdir("data","$(params[:filename]).stl")
  @unpack nmin,nmax,kwargs = params
  args = "\"$filename\",nmin=$nmin,nmax=$nmax,$kwargs"
  Dict(
  "q" => "normal",
  "o" => datadir(jobname*"o.txt"),
  "e" => datadir(jobname*"e.txt"),
  "walltime" => "10:00:00",
  "ncpus" => 1,
  "mem" => "16gb",
  "name" => jobname,
  "func" => params[:func],
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

template = read(scriptsdir("jobtemplate_gadi.sh"),String)

all_params = Dict(
  :func => "rotations_and_displacements",
  :filename => stl_list,
  :nmin => 1,
  :nmax => 14 .* 2 .^ (0:5),
  :kwargs =>  "displacements=exp10.(-17:-1 ), angles=exp10.(-17:-1 )"
  )

for params in dict_list(all_params)
  jobname = savename(params[:filename],params,ignores=(:func,:filename,:kwargs))
  jobname = replace(jobname,"="=>"")
  jobfile = datadir(jobname*".sh")
  open(jobfile,"w") do io
    render(io,template,jobdict(jobname,params))
  end
end





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

include(testdir("data","thingi10k_quality_filter.jl")) # file_ids = Int[]


using STLCutters.Tests: download_things
# avoid re-download
#  save filename, id (check before each download)
#  download by chanks (100 eg) and write jobfile after download from#_to#_
#  save indexes and paths in .bson files
#


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




#@time download_things(file_ids[1:100],path=tmpdir())

## TODO:
#
#  Add HPC name to the ouput
#  Classify folders by HPC
#  Add jobs running all the 5k stls in chanks of 100

