gb2str(gb) = "$(gb)gb"
gb2mb(gb) = gb*1024

gadi_dict = Dict{Symbol,Any}(
  :hpcname => "gadi",
  :gb2mem => gb2str,
  :queue => "normal",
  :julia => "julia" )
  
titani_dict = Dict{Symbol,Any}(
  :hpcname => "titani",                               
  :gb2mem => gb2mb,
  :queue => "serial",
  :julia => "\$HOME/julia-1.5.3/bin/julia" )

acuario_dict = Dict{Symbol,Any}(
  :hpcname => "acuario",                               
  :gb2mem => gb2mb,
  :queue => "R630",
  :julia => "\$HOME/julia-1.5.3/bin/julia" )

hpc_dict = Dict{Symbol,Dict}(
  :gadi => gadi_dict,
  :titani => titani_dict,
  :acuario => acuario_dict )

function jobdict(hpc_id,jobparams,testparams)
  @unpack hpcname,queue,gb2mem,julia= hpc_dict[hpc_id]
  @unpack jobname,walltime,memory,ncpus = jobparams
  @unpack incl,func,args,kwargs = testparams
  Dict(
  "q" => queue,
  "o" => datadir(hpcname,jobname*".out"),
  "e" => datadir(hpcname,jobname*".err"),
  "walltime" => walltime,
  "ncpus" => ncpus,
  "mem" => gb2mem(memory),
  "name" => jobname,
  "julia" => julia,
  "includes" => incl,
  "func" => func,
  "args" => args*kwargs,
  "projectdir" => projectdir()
  )
end


