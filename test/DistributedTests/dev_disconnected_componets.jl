

using STLCutters
using Gridap
using Gridap.Geometry
using Gridap.Arrays
using Gridap.ReferenceFEs
using GridapDistributed
using PartitionedArrays
using GridapEmbedded



function unique_reffe(grid::UnstructuredGrid)
  coords = get_node_coordinates(grid)
  conn = get_cell_node_ids(grid)
  reffes = get_reffes(grid)
  @assert all(==(reffes[1]),reffes)
  ctypes = ones(Int8,length(conn))
  reffes = [reffes[1]]
  UnstructuredGrid(coords,conn,reffes,ctypes)
end

function unique_reffe(model::UnstructuredDiscreteModel)
  grid = get_grid(model)
  grid = unique_reffe(grid)
  UnstructuredDiscreteModel(grid)
end

function unique_reffe(model::GridapDistributed.DistributedDiscreteModel)
  models = map(local_views(model)) do m
    unique_reffe(m)
  end
  GridapDistributed.DistributedDiscreteModel(models,get_cell_gids(model))
end

function generate_components(np,nlc,d=2)
  N = prod(np)
  @assert N % nlc == 0
  comp_to_part = map(1:N) do i
    ( (i-1) ÷ (nlc*d) )*d + (i-1) % d + 1
  end

  part_to_comp = map(1:N) do p
    if p ≤ N ÷ nlc
      findall(==(p),comp_to_part)
    else
      Int[]
    end
  end

  part_to_comp
end

function disconnected_cartesian_model(distribute,components,np,args...)
  # Pseudo distributed model
  ranks = LinearIndices((prod(np),))
  model = CartesianDiscreteModel(ranks,np,args...)
  gids = get_cell_gids(model)
  # Setup exchange components
  graph = reverse(ExchangeGraph(components))
  p = gather(graph.snd,destination=:all) |> PartitionedArrays.getany
  part_to_part = p.data
  lmodels = map(m->[m],local_views(model))
  l_to_g = map(ids->[local_to_global(ids)], partition(gids))
  l_to_p = map(partition(gids)) do ids
    [ map(Reindex(part_to_part),local_to_owner(ids)) ]
  end
  # Send components
  lmodels = exchange(lmodels,graph) |> fetch
  l_to_g = exchange(l_to_g,graph) |> fetch
  l_to_p = exchange(l_to_p,graph) |> fetch
  # Merge components
  models = map(lmodels) do lmodels
    if !isempty(lmodels)
      m = lmodels[1]
      for m_i in lmodels[2:end]
        m = lazy_append(m,m_i) |> UnstructuredDiscreteModel
      end
      m
    end
  end
  l_to_g = map(l_to_g) do l_to_g
    if !isempty(l_to_g)
      _l_to_g = collect(l_to_g[1])
      for l_to_g_i in l_to_g[2:end]
        append!(_l_to_g,l_to_g_i)
      end
      _l_to_g
    else
      eltype(eltype(l_to_g))[]
    end
  end
  l_to_p = map(l_to_p) do l_to_p
    if !isempty(l_to_p)
      _l_to_p = collect(l_to_p[1])
      for l_to_p_i in l_to_p[2:end]
        append!(_l_to_p,l_to_p_i)
      end
      _l_to_p
    else
      eltype(eltype(l_to_p))[]
    end
  end
  lids = map(partition(gids),l_to_g,l_to_p) do ids,l_to_g,l_to_p
    l_to_p = collect(Int32,l_to_p)
    LocalIndices(length(gids),part_id(ids),l_to_g,l_to_p)
  end
  # Create true distributed model
  comp_map = findall(!isempty,components)
  Np = length(comp_map)
  _ranks = distribute(LinearIndices((Np,)))
  _models = map(_ranks) do p
    models[comp_map[p]]
  end
  _lids = map(_ranks) do p
    lids[comp_map[p]]
  end
  _gids = PRange(_lids)
  GridapDistributed.DistributedDiscreteModel(_models,_gids) |> GridapEmbedded.Distributed.merge_nodes |> unique_reffe
end

np = (4,4,4)
nc = (16,16,16)
# domain = (0, 1, 0, 1)
# components = [[1,3],[2,4],[],[]]
distribute = Vector

ranks = distribute(LinearIndices((prod(np),)))
# components = generate_components(np,2,2)
# model = disconnected_cartesian_model(distribute,components,np,domain,nc)


geo = STLGeometry("test/data/cube.stl")

pmin,pmax = get_bounding_box(geo)
pmin -= (pmax-pmin)*0.1
pmax += (pmax-pmin)*0.1

components = generate_components(np,2,2)
model = disconnected_cartesian_model(distribute,components,np,pmin,pmax,nc)
# model = CartesianDiscreteModel(ranks,np,pmin,pmax,nc)


writevtk(Triangulation(model),"disc_model")
writevtk(geo,"geo")



cmodel = CartesianDiscreteModel(pmin,pmax,nc)



cutgeo = cut(model,geo)

Ω = Triangulation(cutgeo,IN)


writevtk(Ω,"trian")

# TODO:
# Implement disconnected_components -> local, global components
# Implement neighbors
# Vectorize propagation
