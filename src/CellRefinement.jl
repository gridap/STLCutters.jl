
const TOL = 1e6*eps()

function insert_vertices!(T,X,p,stl::DiscreteModel,V,f,Tnew,fnew)
  for (k,Vk,fk) in zip(T,V,f)

    v_i = UNSET
    while length(Vk) > 0
      i = farthest_vertex_from_boundary(k,X,p,stl,Vk)
      point = get_vertex(stl,Vk[i])
      dist = distance_to_boundary(k,X,p,point)
      if dist < TOL
        point = move_vertex_to_cell_boundary(k,X,p,point)
        get_vertex_coordinates(get_grid_topology(stl))[Vk[i]] = point
        deleteat!(Vk,i)
      else
        v_i = i
        break
      end
    end

    if length(Vk) > 0 
      v = Vk[v_i]
      deleteat!(Vk,v_i)
      Tk,Xnew = vertex_refinement(k,X,p,get_vertex(stl,v))
      append!(X,Xnew)
      VTk = distribute_vertices(Tk,X,p,stl,Vk)
      fTk = fill([fk;v],length(Tk))
      insert_vertices!(Tk,X,p,stl,VTk,fTk,Tnew,fnew)
    else
      push!(Tnew,k)
      push!(fnew,fk)
    end
  end
end

function insert_edges!(T,X,p,stl::DiscreteModel,E,f,Tnew,fnew,vs)
  for (k,Ek,fk) in zip(T,E,f)

    while length(Ek) > 0
      e = Ek[1]
      edge = get_edge(stl,e)
      if is_on_boundary(k,X,p,edge,atol=TOL)
        popfirst!(Ek)
      else
        break
      end
    end

    if length(Ek) > 0
      e = popfirst!(Ek)
      Tk,Xnew = edge_refinement(k,X,p,get_edge(stl,e),vs)
      append!(X,Xnew)
      ETk = distribute_edges(Tk,X,p,stl,Ek)
      _f =  e+get_offset(get_grid_topology(stl),1)
      fTk = fill([fk;_f],length(Tk))
      insert_edges!(Tk,X,p,stl,ETk,fTk,Tnew,fnew,vs)
    else
      push!(Tnew,k)
      push!(fnew,fk)
    end
  end
end

function insert_facets!(T,X,p,stl,F,f,Tnew,fnew,cell_types,cell_to_io)
  for (k,Fk,fk) in zip(T,F,f)
    if length(Fk) > 0
      Tk,Xnew,c_io,cts = facet_refinement(k,X,p,stl,Fk)
      append!(Tnew,Tk)
      append!(X,Xnew)
      append!(cell_to_io,c_io)
      append!(cell_types,cts)
      _fk = [fk;Fk.+get_offset(get_grid_topology(stl),num_dims(stl))]
      append!(fnew,fill(_fk,length(Tk)))
    else
      push!(Tnew,k)
      push!(cell_to_io,UNSET)
      push!(cell_types,HEX_AXIS)
      push!(fnew,fk)
    end
  end
end

function define_cells!(grid::Grid,f,cell_to_io)
  for i in 1:num_cells(grid)
    if cell_to_io[i] == UNSET
      define_cell!(grid,i,f,cell_to_io)
    end
  end
end

## Helpers

function define_cell!(grid::Grid,i::Integer,f,cell_to_io)
  neig = first_sibling(grid,i,f,cell_to_io)
  if cell_to_io[neig] == UNSET
    define_cell!(grid,neig,f,cell_to_io)
  end 
  cell_to_io[i] = cell_to_io[neig]
end

function first_sibling(grid::Grid,i::Integer,f,c_to_io)
  fi = f[i]
  for k in 1:num_cells(grid)
    if k ≠ i
      fk = f[k]
      if fi ⊆ fk && 
         c_to_io ≠ CellMeshes.IGNORE_FACE &&
         are_connected(grid,i,k,0) && 
         are_in_touch(grid,i,k)

        return k
      end
    end
  end
  @unreachable
end

function are_connected(grid::Grid,i::Integer,k::Integer,d::Integer)
  p_i = get_polytope( get_reffes(grid)[ get_cell_type(grid)[i] ] )
  p_k = get_polytope( get_reffes(grid)[ get_cell_type(grid)[k] ] )
  nodes_i = get_cell_nodes(grid)[i]
  nodes_k = get_cell_nodes(grid)[k]
  faces_i = get_face_vertices(p_i,d)
  faces_k = get_face_vertices(p_k,d)
  for fi in faces_i, fk in faces_k
    connected = true
    for ni in fi, nk in fk
      if nodes_i[ni] ≠ nodes_k[nk]
        connected = false
      end
    end
    if connected
      return true
    end
  end
  false
end

function are_in_touch(grid::Grid,i::Integer,k::Integer)
  cell_i = get_cell(grid,i)
  cell_k = get_cell(grid,k)
  for f in 1:num_facets(cell_k)
    facet = get_facet(cell_k,f)
    if measure(facet) > 0 && is_on_boundary(cell_i,facet,atol=TOL)
      return true
    end
  end
  false
end

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
    if dist ≥ max_dist
      max_dist = dist
      v_i = i
    end
  end
  @assert v_i ≠ UNSET
  v_i
end

function distance_to_boundary(
  cell_nodes,
  node_to_coordinates::Vector{<:Point},
  p::Polytope,
  point::Point)

  pmin,pmax = get_bounding_box(cell_nodes,node_to_coordinates,p)
  abs( min( minimum(point-pmin), minimum(pmax-point) ) )
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
