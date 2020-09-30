
const TOL = 1e6*eps()

function insert_vertices!(T,X,p,stl::DiscreteModel,V,Tnew,Tnew_to_v,v_in)
  for (k,Vk) in zip(T,V)

    v_i = UNSET
    while length(Vk) > 0
      i = farthest_vertex_from_boundary(k,X,p,stl,Vk)
      point = get_vertex(stl,Vk[i])
      dist = distance_to_boundary(k,X,p,point)
      if dist < TOL
        point = move_vertex_to_cell_boundary(k,X,p,point)
        get_vertex_coordinates(get_grid_topology(stl))[i] = point
        deleteat!(Vk,i)
      else
        v_i = i
        break
      end
    end

    if length(Vk) > 0 
      v = Vk[v_i]
      deleteat!(Vk,v_i)
      v_in_k = [ v_in ; [v] ]
      Tk,Xnew = vertex_refinement(k,X,p,get_vertex(stl,v))
      append!(X,Xnew)
      VTk = distribute_vertices(Tk,X,p,stl,Vk)
      insert_vertices!(Tk,X,p,stl,VTk,Tnew,Tnew_to_v,v_in_k)
    else
      push!(Tnew,k)
      push!(Tnew_to_v,v_in)
    end
  end
end

function insert_edges!(T,X,p,stl::DiscreteModel,E,Tnew,Tnew_to_e,e_in,vs)
  for (k,Ek) in zip(T,E)

    if length(Ek) > 0
      e = popfirst!(Ek)
      e_in_k = [ e_in ; [e] ]
      Tk,Xnew = edge_refinement(k,X,p,get_edge(stl,e),vs)
      append!(X,Xnew)
      ETk = distribute_edges(Tk,X,p,stl,Ek)
      insert_edges!(Tk,X,p,stl,ETk,Tnew,Tnew_to_e,e_in_k,vs)
    else
      push!(Tnew,k)
      push!(Tnew_to_e,e_in)
    end
  end
end

function insert_facets!(T,X,p,stl,F,Tnew,cell_types,cell_to_io)
  for (k,Fk) in zip(T,F)
    if length(Fk) > 0
      Tk,Xnew,c_io = facet_refinement(k,X,p,stl,Fk)
      append!(Tnew,Tk)
      append!(X,Xnew)
      append!(cell_to_io,c_io)
      append!(cell_types,fill(TET_AXIS,length(Tk)))
    else
      push!(Tnew,k)
      push!(cell_to_io,UNSET)
      push!(cell_types,HEX_AXIS)
    end
  end
end

## Helpers

function initial_mesh(p::Polytope)
  X = collect(get_vertex_coordinates(p))
  K = collect(1:length(X))
  T = [K]
  T,X
end

function move_vertex_to_cell_boundary(
  cell_nodes::Vector{<:Integer},
  node_to_coordinates::Vector{<:Point},
  p::Polytope,
  point::Point{D}) where D

  pmin,pmax = get_bounding_box(cell_nodes,node_to_coordinates,p)
  for d in 1:D
    if point[d] - pmin[d] < TOL 
      point = Base.setindex(point,pmin[d],d)
    elseif pmax[d] - point[d] < TOL
      point = Base.setindex(point,pmax[d],d)
    end
  end
  point
end

function farthest_vertex_from_boundary(
  cell_nodes::Vector{<:Integer},
  node_to_coordinates::Vector{<:Point},
  p::Polytope,
  stl::DiscreteModel,
  vertices::Vector)

  v_i = UNSET
  max_dist = 0.0
  for (i,v) in enumerate(vertices)
    point = get_vertex(stl,v)
    dist = distance_to_boundary(cell_nodes,node_to_coordinates,p,point)
    if dist > max_dist
      max_dist = dist
      v_i = i
    end
  end
  v_i
end

function distance_to_boundary(cell_nodes,node_to_coordinates,p::Polytope,point::Point)
  @assert have_intersection(cell_nodes,node_to_coordinates,p,point)
  pmin,pmax = get_bounding_box(cell_nodes,node_to_coordinates,p)
  min( minimum(point-pmin), minimum(pmax-point) )
end

function farthest_axis_from_boundary(cell_nodes,node_to_coordinates,p::Polytope,point::Point)
  @assert have_intersection(cell_nodes,node_to_coordinates,p,point)
  pmin,pmax = get_bounding_box(cell_nodes,node_to_coordinates,p)
  max_dists = max( point-pmin, pmax-point )
  _,d = findmax(max_dists.data)
  d
end

function compute_grid(
  cell_to_nodes::AbstractArray,
  node_to_coordinates::Vector{<:Point},
  p::Polytope)

  T = Table(cell_to_nodes)
  X = node_to_coordinates
  reffes = [LagrangianRefFE(p)]
  cell_types = fill(1,length(T))
  UnstructuredGrid(X,T,reffes,cell_types)
end
