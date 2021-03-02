using DrWatson

@quickactivate "STLCutters"

using DataFrames

using STLCutters
using STLCutters.Tests
using STLCutters.Tests: scriptsdir

using Plots

function save_plot(name,data,xfield,yfield;relative=false,labels=Dict())

  if relative
    r0 = filter( i -> i[xfield] == 0, data )
    min_n = minimum(r0.nmax)
    r0 = filter( i -> i.nmax == min_n, r0 )
    y0 = only(r0[!,yfield])
  end
  !isempty(filter( i -> i[xfield] ≠ 0, data )) || return
  
  nmaxs = sort(unique(data.nmax))
  plot()
  for (i,nmax) in enumerate(nmaxs)
    r = filter(i->i.nmax == nmax,data)
    r = filter( i -> i[xfield] ≠ 0, r )
    !isempty(r) || continue
    x = r[!,xfield]
    y = r[!,yfield]
    !relative || ( y = @. abs( y - y0) )
    scatter!(x,y,label="N_max = $nmax",markershape=:auto)
  end
  plot!(xscale=:log10)
  xlabel = haskey(labels,xfield) ? labels[xfield] : "$xfield"
  ylabel = haskey(labels,yfield) ? labels[yfield] : "$yfield"
  plot!(;xlabel,ylabel)

  savefig(plotsdir("$name.pdf"))
end

raw = collect_results(datadir("gadi"))

xfields = [:rotation,:displacement]
yfields = [:surface_error,:volume_error,:min_subcells_x_cell,:max_subcells_x_cell,:avg_subcells_x_cell,:time] 
relative_yfields = [:domain_volume,:domain_surface]

labels = Dict(
  :rotation => "Rotation angle (θ) [rads]", 
  :displacement => "Displacement (Δx)",
  :volume_error => "Domain volume error (ϵ_V)",   
  :domain_volume => "Domain volume variation (Δv_Ω)",
  :surface_error => "Domain surface error (ϵ_Γ)",   
  :domain_surface => "Domain surface variation (ΔΓ_Ω)",
  :time => "Trianglulation time (t) [secs]",
  :num_cells => "Num background cells",
  :num_stl_facets => "Num STL facets" )

all_params = Dict(
  :name => unique(raw[!,:name]),
  :x => xfields,
  :y => [yfields;relative_yfields]  )

xfilter = Dict( :rotation => :displacement, :displacement => :rotation )

for params in dict_list(all_params)
  @unpack name,x,y = params

  r = raw
  r = filter(i->i.name == name,r)
  r = filter(i->i[xfilter[x]] == 0,r)

  plotname = savename(params)
  plotname = replace(plotname,'='=>'_')

  save_plot(plotname,r,x,y,relative=y∈relative_yfields;labels)

end

raw = collect_results(datadir("titani"))

function save_plot(data,xfield,yfield)
  markerkw = (markersize=2,markercolor=:black)

  x = raw[!,xfield]
  y = raw[!,yfield]
  y = @. abs(y) + iszero(y)*eps()/2
  
  scatter(x,y;markerkw...)
  plot!(legend=:none,xscale=:log10,yscale=:log10)
  plot!(xlabel=labels[xfield])
  plot!(ylabel=labels[yfield])
  
  savefig(plotsdir("$(xfield)_$(yfield).pdf"))
end

xfield = :num_stl_facets
yfield = :surface_error
save_plot(raw,xfield,yfield)


xfield = :num_stl_facets
yfield = :volume_error
save_plot(raw,xfield,yfield)
