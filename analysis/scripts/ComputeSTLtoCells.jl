module ComputeSTLtoCells

using DrWatson

@quickactivate "STLCuttersAnalysis"

using Gridap
using STLCutters
using STLCutters: STL
using STLCutters: compute_cell_to_facets
using STLCuttersAnalysis: compute_sizes

function stl_to_cell_stats(stl,bgmesh)
  c_to_stlf = compute_cell_to_facets(bgmesh,stl)
  min_stlf = typemax(Int)
  max_stlf = 0
  sum_stlf = 0
  num_cells_cut = 0
  for cell in 1:length(c_to_stlf)
    nstlf = length(c_to_stlf[cell])
    nstlf > 0 || continue
    max_stlf = max( max_stlf, nstlf )
    min_stlf = min( min_stlf, nstlf )
    sum_stlf += nstlf
    num_cells_cut += 1
  end
  avg_stlf = sum_stlf / num_cells_cut
  min_stlf, max_stlf, avg_stlf, num_cells_cut
end

function get_stl(filename)
  X,T,N = read_stl(filename)
  stl = compute_stl_model(T,X)
  stl = merge_and_collapse(stl)
  writevtk(stl,"stl")
  STL(stl)
end

function get_bgmesh(stl;nmin,nmax,δ)
  pmin,pmax = get_bounding_box(stl)
  Δ = (pmax-pmin)*δ
  origin,sizes,partition = compute_sizes(pmin-Δ,pmax+Δ;nmin,nmax)
  CartesianGrid(origin,sizes,partition)
end

function save_stats(filename;params...)
  name = first(splitext(basename(filename)))
  title = savename(name,params)
  println("Running $title ...")

  stl = get_stl(filename)
  bgmesh = get_bgmesh(stl;params...)
  t = @timed min,max,avg,ncut = stl_to_cell_stats(stl,bgmesh)

  min_stl_faces_x_cell = min
  max_stl_faces_x_cell = max
  avg_stl_faces_x_cell = avg
  num_cells_cut = ncut

  println("TIME = $(t.time) s")
  println("Num Cut Cells : $ncut")
  println("Min Num STL Faces x Cell : $min")
  println("Max Num STL Faces x Cell : $max")
  println("Avg Num STL Faces x Cell : $avg")

  out = Dict{Symbol,Any}()
  @pack! out = name, num_cells_cut
  @pack! out = min_stl_faces_x_cell, max_stl_faces_x_cell, avg_stl_faces_x_cell
  push!(out,params...)

  outname = datadir("local",title*".bson")
  save(outname,out)
end

stl_list = [
  "cube",
  "37881", # Gear
  "47076", # Azteca pyramid
  "96457", # NightHawk
  "550964", # Eiffel Tower
  "293137", # Bunny low-poly
  "441708", # Standford Bunny
  "35269", # Octocat
  "65904_1", # Heart 
  "551021", # Arc de Triomph
  "37266", # Extruded Earth
  "252119"] # Angel

nmaxs = 14 * 2 .^ (0:5)

all_params = Dict(
  :name => stl_list,
  :nmax => nmaxs,
  :nmin => 1,
  :δ => 0.2
 )

for params in dict_list(all_params)
  @unpack name,nmax,nmin,δ = params
  filename = datadir("geometries",name*".stl")
  save_stats(filename;nmax,nmin,δ)
end




end # modul
