gb2str(gb) = "$(gb)gb"
gb2mb(gb) = gb*1024



gadi_dict = Dict{Symbol,Any}(
  :hpcname => "gadi",
  :memory => gb2str,
  :queue => "normal",
  :julia => "julia" )
  
titani_dict = Dict{Symbol,Any}(
  :hpcname => "titani",                               
  :memory => gb2mb,
  :queue => "serial",
  :julia => "\$HOME/julia-1.5.3/bin/julia" )

acuario_dict = Dict{Symbol,Any}(
  :hpcname => "acuario",                               
  :memory => gb2mb,
  :queue => "R630",
  :julia => "\$HOME/julia-1.5.3/bin/julia" )

hpc_dict = Dict{Symbol,Dict}(
  :gadi => gadi_dict,
  :titani => titani_dict,
  :acuario => acuario_dict )
