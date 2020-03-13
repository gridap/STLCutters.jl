
struct BulkMesh{D,T,M}
  background_mesh::M
  cell_to_vertex_coordinates::Table{Point{D,T}}
  cell_to_d_to_dface_to_vertices::Table{Table{Int}}
  cell_to_d_to_dface_to_in_out_boundary::Vector{Table{Int}}
  cell_to_in_out_cut::Vector{Int}
  vertex_to_in_out_boundary::Vector{Int}

  # TODO: add constructor -> (bg_mesh,sm)
end

num_dims(m::BulkMesh{D}) where D = D

get_background_mesh(m::BulkMesh) = m.background_mesh

num_cells(m::BulkMesh) = num_cells( get_background_mesh(m) )

num_vertices(m::BulkMesh) = num_vertices( get_background_mesh(m) )

num_vertices_per_cell(m::BulkMesh) = num_vertices_per_cell( get_background_mesh(m) )

function get_vertex_id(m::BulkMesh,cell::Integer,lvertex::Integer)
  bg_mesh = get_background_mesh(m)
  get_vertex_id(bg_mesh,cell,lvertex)
end

get_cell_in_out(m::BulkMesh,cell::Integer) = m.cell_to_in_out_cut[cell]

is_cell_cut(m::BulkMesh,cell::Integer) = m.cell_to_in_out_cut[cell] == FACE_CUT

function get_vertex_coordinates(m::BulkMesh,i::Integer)
  bg_mesh = get_background_mesh(m)
  get_vertex_coordinates(bg_mesh,i)
end

function get_dface_to_vertices(m::BulkMesh,cell::Integer,d::Integer)
  m.cell_to_d_to_dface_to_vertices[cell,d+1]
end

function num_subfaces(m::BulkMesh,cell::Integer,d::Integer)
  df_to_v = get_dface_to_vertices(m,cell,d)
  length(df_to_v)
end

function num_subfaces(m::BulkMesh)
  n = 0
  for d in 0:num_dims(m)
    n += num_subfaces(m,d)
  end
  n
end

function num_subfaces(m::BulkMesh,d::Integer)
  n = 0
  for k in 1:num_cells(m)
    if is_cell_cut(m,k)
      n += num_subfaces(m,k,d)
    end
  end
  n
end

num_subvertices(m::BulkMesh) = num_subfaces(m,0)

num_subvertices(m::BulkMesh,cell::Integer) = num_subfaces(m,cell,0)

num_subedges(m::BulkMesh) = num_subfaces(m,1)

num_subedges(m::BulkMesh,cell::Integer) = num_subfaces(m,cell,1)

num_subfacets(m::BulkMesh{D}) where D = num_subfaces(m,D-1)

num_subfacets(m::BulkMesh,cell::Integer) where D = num_subfaces(m,cell,D-1)

num_subcells(m::BulkMesh{D}) where D = num_subfaces(m,D)

num_subcells(m::BulkMesh{D},cell::Integer) where D = num_subfaces(m,cell,D)

function get_subvertex_coordinates(m::BulkMesh,cell::Integer,vertex::Integer)
  m.cell_to_vertex_coordinates[cell,vertex]
end

function get_subface_in_out_boundary(m::BulkMesh,cell::Integer,d::Integer,face::Integer)
  m.cell_to_d_to_dface_to_in_out_boundary[cell][d+1,face]
end

function writevtk(m::BulkMesh{D,T},file_base_name) where {D,T}
  d_to_vtk_tet_type_id = Dict(0=>1,1=>3,2=>5,3=>10)
  d_to_vtk_hex_type_id = Dict(0=>1,1=>3,2=>9,3=>12)
  vtk_hex_lvertices = [1,2,4,3,5,6,8,7]

  points = zeros(T,D,num_subvertices(m)+num_vertices(m))
  for i in 1:num_vertices(m)
    p = get_vertex_coordinates(m,i)
    for d in 1:D
      points[d,i] = p[d]
    end
  end
  offset = num_vertices(m)
  for k in 1:num_cells(m)
    if is_cell_cut(m,k)
      for i in 1:num_subvertices(m,k)
        p = get_subvertex_coordinates(m,k,i)
        for d in 1:D
          points[d,offset+i] = p[d]
        end
      end
      offset += num_subvertices(m,k)
    end
  end
  
  cells = MeshCell{Vector{Int64}}[]
  offset = num_vertices(m)
  for k in 1:num_cells(m)
    if is_cell_cut(m,k)
      for d in D:D
        dface_to_vertices = get_dface_to_vertices(m,k,d)
        vtk_type = VTKCellType(d_to_vtk_tet_type_id[d])
        for i in 1:length(dface_to_vertices)
          vertices = [ offset+dface_to_vertices[i,j] for j in 1:length(dface_to_vertices,i) ]
          push!( cells, MeshCell(vtk_type,vertices) )
        end
      end
      offset += num_subvertices(m,k)
    else
      vtk_type = VTKCellType(d_to_vtk_hex_type_id[D])
      vertices = zeros(Int,num_vertices_per_cell(m))
      for j in 1:num_vertices_per_cell(m)
        lv = vtk_hex_lvertices[j]
        vertices[j] = get_vertex_id(m,k,lv)
      end
      push!( cells, MeshCell(vtk_type,vertices) )
    end
  end

  vtkfile = vtk_grid(file_base_name,points,cells)
  vtkfile["IO",VTKCellData()] = _get_in_out_face_data(m,D:D)
  vtk_save(vtkfile)
end

function _get_in_out_face_data(m::BulkMesh,range)
  data = Int[]
  for k in 1:num_cells(m)
    if is_cell_cut(m,k)
      for d in range, i in 1:num_subfaces(m,k,d)
        push!(data,get_subface_in_out_boundary(m,k,d,i))
      end
    else
      push!(data,get_cell_in_out(m,k))
    end
  end
  data
end



