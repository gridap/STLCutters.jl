module PlotTestMatrix

using Plots
using CSV
using DataFrames

data = CSV.read("out_test_matrix.csv",DataFrame)

geometries = unique(data.geometry)


for geom in geometries
  data_g = filter( x -> x.geometry == geom, data )
  data_g = filter( x -> x.kdtree == false, data )
  N = unique(data_g.n)
  plot(markershape=:auto)
  for n in N
    data_n = filter( x -> x.n == n, data_g )
    data_n = filter( x -> x.displacement == 0 && x.rotation > 0, data_n )
    errors = data_n.domain_error
    scatter!(data_n.rotation,abs.(errors[1].-errors),label="n=$n")
  end
  plot!(xscale=:log10)
  plot!(xlabel="Rotation angle (θ) [rads]")
  plot!(ylabel="Domain volume variation (Δv_Ω)")
  savefig("$(geom)_domain_error_rotation.pdf")
end

for geom in geometries
  data_g = filter( x -> x.geometry == geom, data )
  data_g = filter( x -> x.kdtree == false, data )
  N = unique(data_g.n)
  plot(markershape=:auto)
  for n in N
    data_n = filter( x -> x.n == n, data_g )
    data_n = filter( x -> x.displacement == 0 && x.rotation > 0, data_n )
    scatter!(data_n.rotation,data_n.time,label="n=$n")
  end
  plot!(xscale=:log10)
  plot!(xlabel="Rotation angle (θ) [rads]")
  plot!(ylabel="Triangulation time (t) [secs]")
  savefig("$(geom)_triangulation_time.pdf")
end

for geom in geometries
  data_g = filter( x -> x.geometry == geom, data )
  data_g = filter( x -> x.displacement == 0 && x.rotation == 0, data_g )
  t = plot(markershape=:auto)
  s = plot(markershape=:auto)
  m = plot(markershape=:auto)
  for kdtree in (true,false)
    data_k = filter( x -> x.kdtree == kdtree, data_g )
    scatter!(t,data_k.n,data_k.time,label="Kd-Tree ($kdtree)")
    scatter!(s,data_k.n,data_k.num_subcells,label="Kd-Tree ($kdtree)")
    scatter!(m,data_k.n,data_k.max_num_subcells,label="Kd-Tree ($kdtree)")
  end
  plot!(t,xlabel="Mesh size (n)")
  plot!(t,ylabel="Triangulation time (t) [secs]")
  savefig(t,"$(geom)_kd_tree_effect_time.pdf")
  plot!(s,xlabel="Mesh size (n)")
  plot!(s,ylabel="Num subcells (Ns)")
  savefig(s,"$(geom)_kd_tree_effect_subcells.pdf")
  plot!(m,xlabel="Mesh size (n)")
  plot!(m,ylabel="Max num subcells ( max(Ns) )")
  savefig(m,"$(geom)_kd_tree_effect_max_subcells.pdf")
end

end # module
