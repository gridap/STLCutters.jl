using DrWatson

@quickactivate "STLCuttersAnalisis"

using DataFrames
using Plots

# Functions

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
    if xfield != :nmax || yfield ∉ (:time,:error_l2,:error_h1) 
      x = x[2:end]
      y = y[2:end]
    end
    if  xfield == :nmax && yfield ∈ (:error_h1,:error_l2)
      x .*= 2
    end
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
  (xfield == :nmax) && plot!(xticks=( 2 .^ (1:5), 2 .^ (1:5) ))
  plot!(xlabel=labels[xfield],ylabel=labels[yfield])
  (yfield == :time) && plot!(ylabel="Relative time (t/t₀)")
  if xfield == :nmax && yfield ∈ (:error_h1,:error_l2)
    xlims = Plots.get_sp_lims(plot!()[1],:x)
    ylims = Plots.get_sp_lims(plot!()[1],:y)
    order = yfield == :error_h1 ? 1 : 2
    x = [1,1e2]
    y = @. (1/x)^order
    plot!(x,y.*1e-1,color=:black,linestyle=:dash,label="slope = -$order")
    plot!(x,y.*min_y0 .* 2. ^ (order-2),color=:black,linestyle=:dash,label=:none)
    plot!(;xlims,ylims)
  end
  savefig(plotsdir("$plotname.pdf"))
end

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

# Filter Raw Data

raw_gadi = collect_results(datadir("gadi"))
raw_titani = collect_results(datadir("titani"))

raw_strong_scaling = filter(i->!ismissing(i.nruns)&&i.nruns>1,raw_gadi)

raw_gadi = filter(i->ismissing(i.nruns)||i.nruns==1,raw_gadi)

raw_matrix = filter(i->i.nmax ≠ 100,raw_gadi)

raw_10k_gadi = filter(i->i.nmax == 100,raw_gadi)
raw_10k_titani = raw_titani

raw_10k = vcat( raw_10k_gadi,raw_10k_titani,cols=:union)

raw_10k = filter( i->!contains(i.path,"#"), raw_10k)

raw_matrix = filter(i->i.name≠"96457"||i.rotation<0.1,raw_matrix)
replace!(raw_matrix.name,"65904_1"=>"65904")
filter!(i->i.name≠"wine_glass",raw_matrix)

raw_matrix_poisson = filter( i-> !ismissing(i.poisson) && i.poisson, raw_matrix)
raw_matrix = filter( i -> ismissing(i.poisson) || !i.poisson, raw_matrix)

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


# Rotations and displacements

test_fields = [:rotation,:displacement,:nmax]
default_value = Dict(
  :rotation => 0,
  :displacement => 0,
  :nmax => 112 )
data = raw_matrix

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

# Rotations and displacements Poisson

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

# Error norm convergence

data = filter(i->i.solution_order==2,raw_matrix_poisson)
data = filter(i-> !ismissing(i.solver) && i.solver==:amg, data)
data = filter(i-> i.nmax > 14, data)

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

## Thingi10k

include( datadir("geometries","thingi10k_solid_manifold.jl"))

raw_10k = filter(i->parse(Int,i.name)∈file_ids,raw_10k)

# Dot plots 

xfield = :num_stl_facets
yfield = :surface_error
save_plot(raw_10k,xfield,yfield)

xfield = :num_stl_facets
yfield = :volume_error
save_plot(raw_10k,xfield,yfield)

raw_10k_t = filter(i->ismissing(i.solver),raw_10k)
xfield = :num_stl_facets
yfield = :time
save_plot(raw_10k_t,xfield,yfield)

# Histograms

accumulated_frequency_histogram(raw_10k,:volume_error)

accumulated_frequency_histogram(raw_10k,:surface_error)

accumulated_frequency_histogram(raw_10k_t,:time,range=exp10.(0:0.1:6))

accumulated_frequency_histogram(raw_10k,:num_stl_facets,range=exp10.(0:0.1:6))

# Print large errors

tol = 1e-9
n = size(raw_10k,1)
c = count(i->abs(i) < tol,raw_10k.surface_error)
println("$(c/n*100)% of $n has ϵ_Γ below $tol ($c out of $n)")

c = count(i->abs(i) < tol,raw_10k.volume_error)
println("$(c/n*100)% of $n has ϵ_V below $tol ($c out of $n)")

# Strong Scaling

data = raw_strong_scaling

geos = unique(data.name)
nmaxs = sort(unique(data.nmax))
#nthreads = sort(unique(data.nthreads))
#nmaxs = [100]
geos = ["cube"]

plot()
for g in geos, nmax in nmaxs, method in (:threads,:spawn)
  d = filter(i -> i.name == g && i.nmax == nmax && i.threading == method,data)
  nthreads = sort(unique(d.nthreads))
  times = zeros(length(nthreads))
  for (i,n) in enumerate(nthreads)
    _d = filter(i->i.nthreads==n,d)
    min_t = minimum(_d.time)
    times[i] = min_t
  end
  speedup = times[1] ./ times
  plot!(nthreads,speedup,label="n=$nmax,$g,@$method",markershape=:auto)
end

xlims = Plots.get_sp_lims(plot!()[1],:x)
ylims = Plots.get_sp_lims(plot!()[1],:y)
plot!(1:48,1:48,color=:black,linestyle=:dash,label=:none)
plot!(;xlims,ylims)
plot!(legend=:outerbottomright)
plot!(xlabel="Num threads",ylabel="Speed-up")
#savefig(plotsdir("strong_scaling.pdf"))

