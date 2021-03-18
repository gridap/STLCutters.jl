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

raw_10k = vcat( raw_10k_gadi,raw_10k_titani,cols=:union)

raw_matrix = filter(i->i.name≠"96457"||i.rotation<0.1,raw_matrix)

raw_matrix_poisson = filter( i-> !ismissing(i.poisson) && i.poisson, raw_matrix )
raw_matrix = filter( i -> ismissing(i.poisson) || !i.poisson, raw_matrix )

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
  :num_stl_facets => "Num STL facets",
  :nmax => "Relative cell size (h₀/h) ",
  :error_h1 => "H1 error norm (ϵ_H¹)",
  :error_l2 => "L2 error norm (ϵ_L²)" )

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

#  save_plot(plotname,r,x,y;relative,labels,ylog)

end

function plot_all_geometries(raw,plotname,xfield,yfield)
  data = filter(i->!ismissing(i[xfield]),raw)
  for field in test_fields
    field ≠ xfield || continue
    data = filter(i->i[field] == default_value[field],data)
  end
  names = sort(unique(data[!,:name]))
  min_y0 = Inf
  sort!(names)
  plot()
  for name in names
    d = filter(i->i.name == name,data)
    sort!(d,xfield)
    x = d[!,xfield]
    y = d[!,yfield]
    x0,y0 = x[1],y[1]
    min_y0 = min(y0,min_y0)
    (xfield == :nmax ) || ( x = x[2:end] )
    (xfield == :nmax ) || ( y = y[2:end] )
    if yfield ∉ (:error_l2,:error_h1)
      y = yfield == :time ? y / y0 : @. abs( y - y0 ) / y0 
    end
    y = @. y + iszero(y)*eps()/2
    (xfield == :nmax) && ( x = x ./ x0 ) 
    plot!(x,y,label="$name",markershape=:auto)
  end
  plot!(xscale=:log10)
  plot!(yscale=:log10)
  plot!(legend=:outerbottomright)
  (xfield == :nmax) && plot!(xticks=( 2 .^ (0:5), 2 .^ (0:5) ))
  plot!(xlabel=labels[xfield],ylabel=labels[yfield])
  (yfield == :time) && plot!(ylabel="Relative time (t/t₀)")
  if xfield == :nmax && yfield ∈ (:error_h1,:error_l2)
    xlims = Plots.get_sp_lims(plot!()[1],:x)
    ylims = Plots.get_sp_lims(plot!()[1],:y)
    order = yfield == :error_h1 ? 1 : 2
    x = [1,1e2]
    y = @. (1/x)^order
    plot!(x,y.*1e-1,color=:black,linestyle=:dash,label="slope = -$order")
    plot!(x,y.*min_y0./2,color=:black,linestyle=:dash,label=:none)
    plot!(;xlims,ylims)
  end
  savefig(plotsdir("$plotname.pdf"))
end

test_fields = [:rotation,:displacement,:nmax]
default_value = Dict(
  :rotation => 0,
  :displacement => 0,
  :nmax => 112 )

data = raw_matrix
replace!(data.name,"wine_glass"=>"3201401")

all_params = Dict(
  :x => test_fields,
  :y => [:domain_surface,:domain_volume,:time]
)

for params in dict_list(all_params)
  plotname = savename(params)
  plotname = replace(plotname,'='=>'_')

  @unpack x,y = params
  plot_all_geometries(data,plotname,x,y)
end

data = filter(i->i.solution_order==1,raw_matrix_poisson)
test_fields = [:rotation,:displacement]
all_params = Dict(
  :x => test_fields,
  :y => [:error_h1,:error_l2]
)

for params in dict_list(all_params)
  plotname = savename(params)
  plotname = replace(plotname,'='=>'_')

  @unpack x,y = params
  plot_all_geometries(data,plotname,x,y)
end


data = filter(i->i.solution_order==2,raw_matrix_poisson)
test_fields = [:nmax]
all_params = Dict(
  :x => test_fields,
  :y => [:error_h1,:error_l2]
)

for params in dict_list(all_params)
  plotname = savename(params)
  plotname = replace(plotname,'='=>'_')

  @unpack x,y = params
  plot_all_geometries(data,plotname,x,y)
end


include(testdir("data","thingi10k_quality_filter.jl"))

raw_10k = filter(i->parse(Int,i.name)∉blacklist_ids,raw_10k)


function save_plot(data,xfield,yfield;prefix=nothing)
  markerkw = (markersize=2,markercolor=:black)

  x = data[!,xfield]
  y = data[!,yfield]
  y = @. abs(y) + iszero(y)*eps()/2
  
  scatter(x,y;markerkw...)
  plot!(legend=:none,xscale=:log10,yscale=:log10)
  plot!(xlabel=labels[xfield])
  plot!(ylabel=labels[yfield])
  
  filename = "$(xfield)_$(yfield).pdf"
  isnothing(prefix) || ( filename = "$(prefix)_$(filename)" )
  savefig(plotsdir(filename))
end

function accumulated_frequency_histogram(data,field;range=exp10.(-17:0.1:0),prefix=nothing)
  f = @. data[!,field]
  x = []
  counts = []
  for lim in range
    c = count(i->abs(i) > lim,f)
    push!(counts,c)
    push!(x,lim)
    c > 0 || break
  end
  n = length(f)
  y = @. ( 1 - counts / n ) * 100
  plot()
  plot(x,y,linecolor=:black)
  plot!(legend=:none,xscale=:log10)
  plot!(ylims=[0,100])
  plot!(ylabel="Accumulated percentage [%]")
  plot!(xlabel=labels[field])
  filename = "histogram_$field.pdf"
  isnothing(prefix) || ( filename = "$(prefix)_$(filename)" )
  savefig(plotsdir(filename))
end

xfield = :num_stl_facets
yfield = :surface_error
save_plot(raw_10k,xfield,yfield)


xfield = :num_stl_facets
yfield = :volume_error
save_plot(raw_10k,xfield,yfield)

xfield = :num_stl_facets
yfield = :time
save_plot(raw_10k,xfield,yfield)

accumulated_frequency_histogram(raw_10k,:volume_error)

accumulated_frequency_histogram(raw_10k,:surface_error)

accumulated_frequency_histogram(raw_10k,:time,range=exp10.(0:0.1:6))

accumulated_frequency_histogram(raw_10k,:num_stl_facets,range=exp10.(0:0.1:6))


tol = 1e-9
c = count(i->abs(i) < tol,raw_10k.surface_error)
println("$(c/size(raw_10k,1)*100)% of $(size(raw_10k,1)) has ϵ_Γ below $tol ($c out of $(size(raw_10k,1)))")

c = count(i->abs(i) < tol,raw_10k.volume_error)
println("$(c/size(raw_10k,1)*100)% of $(size(raw_10k,1)) has ϵ_V below $tol ($c out of $(size(raw_10k,1)))")

raw_10k = filter(i->ismissing(i.min_h)||i.min_h>1e-5,raw_10k)

xfield = :num_stl_facets
yfield = :surface_error
save_plot(raw_10k,xfield,yfield,prefix="filter")

xfield = :num_stl_facets
yfield = :volume_error
save_plot(raw_10k,xfield,yfield,prefix="filter")

accumulated_frequency_histogram(raw_10k,:volume_error,prefix="filter")

accumulated_frequency_histogram(raw_10k,:surface_error,prefix="filter")

