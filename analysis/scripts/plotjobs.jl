using DrWatson

@quickactivate "STLCutters"

using DataFrames

using STLCutters
using STLCutters.Tests
using STLCutters.Tests: scriptsdir

using Plots

function save_plot(name,data,xfield,yfield;relative=false,labels=Dict(),xlog=true,ylog=false)

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
    !relative || ( y = @. abs( y - y0) / y0 )
    scatter!(x,y,label="N_max = $nmax",markershape=:auto)
  end
  !xlog || plot!(xscale=:log10)
  !ylog || plot!(yscale=:log10)
  plot!(legend=:outerbottomright)
  xlabel = haskey(labels,xfield) ? labels[xfield] : "$xfield"
  ylabel = haskey(labels,yfield) ? labels[yfield] : "$yfield"
  plot!(;xlabel,ylabel)
  plot!(title=name)

  savefig(plotsdir("$name.pdf"))
end

raw_gadi = collect_results(datadir("gadi"))
raw_titani = collect_results(datadir("titani"))

raw_matrix = filter(i->i.nmax ≠ 100,raw_gadi)

raw_10k_gadi = filter(i->i.nmax == 100,raw_gadi)
raw_10k_titani = raw_titani
raw_10k = vcat( raw_10k_gadi,raw_10k_titani)

raw_matrix = filter(i->i.name≠"96457"||i.rotation<0.1,raw_matrix)

xfields = [:rotation,:displacement]
yfields = [:surface_error,:volume_error,:min_subcells_x_cell,:max_subcells_x_cell,:avg_subcells_x_cell,:time] 
relative_yfields = [:domain_volume,:domain_surface]
log_yfields = [:min_subcells_x_cell,:max_subcells_x_cell,:avg_subcells_x_cell,:time]

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
  :name => unique(raw_matrix[!,:name]),
  :x => xfields,
  :y => [yfields;relative_yfields]  )

xfilter = Dict( :rotation => :displacement, :displacement => :rotation )

for params in dict_list(all_params)
  @unpack name,x,y = params

  r = raw_matrix
  r = filter(i->i.name == name,r)
  r = filter(i->i[xfilter[x]] == 0,r)

  plotname = savename(params)
  plotname = replace(plotname,'='=>'_')

  relative = y ∈ relative_yfields
  ylog = y ∈ log_yfields

  save_plot(plotname,r,x,y;relative,labels,ylog)

end

include(testdir("data","thingi10k_quality_filter.jl"))

raw_10k = filter(i->parse(Int,i.name)∉blacklist_ids,raw_10k)
raw_10k = filter(i->ismissing(i.min_h)||i.min_h>1e-5,raw_10k)


function save_plot(data,xfield,yfield)
  markerkw = (markersize=2,markercolor=:black)

  x = data[!,xfield]
  y = data[!,yfield]
  y = @. abs(y) + iszero(y)*eps()/2
  
  scatter(x,y;markerkw...)
  plot!(legend=:none,xscale=:log10,yscale=:log10)
  plot!(xlabel=labels[xfield])
  plot!(ylabel=labels[yfield])
  
  savefig(plotsdir("$(xfield)_$(yfield).pdf"))
end

xfield = :num_stl_facets
yfield = :surface_error
save_plot(raw_10k,xfield,yfield)


xfield = :num_stl_facets
yfield = :volume_error
save_plot(raw_10k,xfield,yfield)
