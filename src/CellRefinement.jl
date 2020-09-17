
const TOL = 1e6*eps()

function insert_vertices!(T,X,p,V,Tnew,STL_vertices,Tnew_to_v,v_in)
  for (k,Vk) in zip(T,V)

    iv = UNSET
    while length(Vk) > 0
      iv = farthest_vertex_from_boundary(k,X,Vk,STL_vertices)
      point = STL_vertices[Vk[iv]]
      dist = distance_to_boundary(k,X,point)
      if dist < TOL
        point = move_vertex_to_cell_boundary(k,X,point)
        STL_vertices[Vk[iv]] = point
        deleteat!(Vk,iv)
      else
        break
      end
    end

    if length(Vk) > 0
      v = Vk[iv]
      deleteat!(Vk,iv)
      v_in_k = [ v_in ; [v] ]
      Tk,Xnew = vertex_refinement(k,X,p,STL_vertices[v])
      append!(X,Xnew)
      VTk = distribute_faces(Tk,X,p,Vk,STL_vertices)
      insert_vertices!(Tk,X,p,VTk,Tnew,STL_vertices,Tnew_to_v,v_in_k)
    else
      push!(Tnew,k)
      push!(Tnew_to_v,v_in)
    end
  end
end

function insert_edges!(T,X,p,E,Tnew,STL_edges,Tnew_to_e,e_in,vs)
  for (k,Ek) in zip(T,E)

    if length(Ek) > 0
      e = popfirst!(Ek)
      e_in_k = [ e_in ; [e] ]
      Tk,Xnew = edge_refinement(k,X,p,STL_edges[e],vs)
      append!(X,Xnew)
      VTk = distribute_faces(Tk,X,p,Ek,STL_edges)
      insert_edges!(Tk,X,p,VTk,Tnew,STL_edges,Tnew_to_e,e_in_k,vs)
    else
      push!(Tnew,k)
      push!(Tnew_to_e,e_in)
    end
  end
end

function insert_facets!(T,X,p,F,Tnew,cell_types,cell_to_io,STL_facets)
  for (k,Fk) in zip(T,F)
    if length(Fk) > 0
      Tk,Xnew,c_io = facet_refinement(k,X,p,Fk,STL_facets)
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
  point::Point{D}) where D

  pmin,pmax = get_bounding_box(cell_nodes,node_to_coordinates)
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
  vertices::Vector,
  STL_vertices)

  iv = UNSET
  max_dist = 0.0
  for (i,v) in enumerate(vertices)
    p = STL_vertices[v]
    dist = distance_to_boundary(cell_nodes,node_to_coordinates,p)
    if dist > max_dist
      max_dist = dist
      iv = i
    end
  end
  iv
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
