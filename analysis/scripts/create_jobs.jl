using Mustache
using DrWatson

function jobdict(jobname,params)
  filename = testdir("data","$(params[:filename]).stl")
  Dict(
  "q" => "normal",
  "o" => datadir(jobname*"o.txt"),
  "e" => datadir(jobname*"e.txt"),
  "walltime" => "01:00:00",
  "ncpus" => 1,
  "mem" => "16gb",
  "name" => jobname,
  "func" => params[:func],
  "filename" => filename,
  "nmin" => params[:nmin],
  "nmax" => params[:nmax],
  "params" => params[:params],
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
  :params =>  ",displacements=exp10.(-17:-1 ), angles=exp10.(-17:-1 )"
  )

for params in dict_list(all_params)
  jobname = savename(params[:filename],params,ignores=(:func,:filename,:params))
  jobname = replace(jobname,"="=>"")
  jobfile = datadir(jobname*".sh")
  open(jobfile,"w") do io
    render(io,template,jobdict(jobname,params))
  end
end

