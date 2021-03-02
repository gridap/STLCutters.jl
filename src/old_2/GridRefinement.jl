
const FACE_IN = CellMeshes.FACE_IN

const FACE_OUT = CellMeshes.FACE_OUT

const FACE_CUT = CellMeshes.FACE_CUT


function refine_grid(grid::Grid,stl::DiscreteModel)
  _,cell_to_stlf = compute_face_to_cells(grid,stl)
  refine_grid(grid,stl,cell_to_stlf)
end

function refine_grid(grid::Grid{D},stl::DiscreteModel,cell_to_stlf) where D

  bgcell_to_ioc = fill( Int8(UNSET), num_cells(grid) )
  bgnode_to_ioc = fill( Int8(UNSET), num_nodes(grid) )
  T = Vector{Int}[]
  X = empty(get_vertex_coordinates(get_grid_topology(stl)))
  if D == 2
    reffes = [QUAD4,TRI3]
  elseif D == 3
    reffes = [HEX8,TET4]
  end
  cell_types = Int8[]
  cell_to_io = Int8[]
  cell_to_bgcell = Int[]
  for cell in 1:num_cells(grid)
    Fs = get_face_lists(stl,cell,cell_to_stlf)
    if D == 2
      V,F = Fs
    elseif D == 3
      V,E,F = Fs
    end
    if length(F) > 0
      Tk,Xk,p = init_cell_mesh(grid,cell) 
      fk = [Int[]]
      Tknew,fknew = empty(Tk), empty(fk)
      VTk = distribute_vertices(Tk,Xk,p,stl,V)
      insert_vertices!(Tk,Xk,p,stl,VTk,fk,Tknew,fknew)
      Tk,fk = Tknew,fknew
      if D == 3
        Tknew,fknew = empty(Tk), empty(fk)
        ETk = distribute_edges(Tk,Xk,p,stl,E)
        vs = get_default_directions(E,stl)
        insert_edges!(Tk,Xk,p,stl,ETk,fk,Tknew,fknew,vs)
        Tk,fk = Tknew,fknew
      end
      Tknew,fknew = empty(Tk), empty(fk)
      k_cell_types,k_cell_to_io = Int8[], Int8[]
      FTk = distribute_facets(Tk,Xk,p,stl,F)
      insert_facets!(Tk,Xk,p,stl,FTk,fk,Tknew,fknew,k_cell_types,k_cell_to_io)
      Tk,fk = Tknew,fknew
      grid_k = UnstructuredGrid(Xk,Table(Tk),reffes,k_cell_types)
      define_cells!(grid_k,fk,k_cell_to_io)
      
      bgcell_to_ioc[cell] = CellMeshes.FACE_CUT
      for f in F
        facet = get_cell(stl,f)
        for bgnode in get_cell_nodes(grid)[cell]
          point = get_node_coordinates(grid)[bgnode]
          if distance(facet,point) < TOL
            bgnode_to_ioc[bgnode] = CellMeshes.FACE_CUT
          end
        end
      end
      for (i,K) in enumerate(Tk), v in K
        if v ≤ num_vertices(p)
          bgnode = get_cell_nodes(grid)[cell][v]
          if bgnode_to_ioc[bgnode] == UNSET
            bgnode_to_ioc[bgnode] = k_cell_to_io[i]
          end
        end
      end
      Tk = [ K.+length(X) for K in Tk ]
      append!(X,Xk)
      append!(T,Tk)
      append!(cell_types,k_cell_types)
      append!(cell_to_io,k_cell_to_io)
      append!(cell_to_bgcell,fill(cell,length(Tk)))
    end
  end

  grid_topology = GridTopology(grid)
  stack = Int[]
  for cell in 1:num_cells(grid)
    if bgcell_to_ioc[cell] == FACE_CUT
      resize!(stack,0)
      push!(stack,cell)
      while length(stack) > 0
        current_cell = pop!(stack)
        for node in get_cell_nodes(grid)[current_cell]
          if bgnode_to_ioc[node] == FACE_IN
            for neig_cell in get_faces(grid_topology,0,num_dims(grid))[node]
              if bgcell_to_ioc[neig_cell] == UNSET
                bgcell_to_ioc[neig_cell] = FACE_IN
                for neig_node in get_cell_nodes(grid)[neig_cell]
                  if bgnode_to_ioc[neig_node] == UNSET
                    bgnode_to_ioc[neig_node] = FACE_IN
                  end
                  push!(stack,neig_cell)
                end
              end
            end
          end
        end
      end
    end
  end
  for cell in 1:num_cells(grid)
    if bgcell_to_ioc[cell] == UNSET
      bgcell_to_ioc[cell] = FACE_OUT
    end
  end
  T,X,reffes,cell_types,cell_to_io,cell_to_bgcell,bgcell_to_ioc
end

function get_face_lists(stl::DiscreteModel{D},cell::Integer,cell_to_stl_faces) where D
  topo = get_grid_topology(stl)
  Fs = ntuple( i -> Int[], Val{D+1}() )
  for face in cell_to_stl_faces[cell]
    d = get_facedims(topo)[face]
    dface = face - get_offset(topo,d)
    push!( Fs[d+1], dface ) 
  end
  Fs
end

function init_cell_mesh(grid::Grid,cell::Integer)
  p = get_polytope(get_cell_reffes(grid)[cell])
  X = get_cell_coordinates(grid)[cell]
  T = [ collect( 1:num_vertices(p) ) ]
  T,X,p
end

