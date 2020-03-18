
struct SubTriangulation{D,T}
  cell_to_points::Table{Int}
  cell_to_inoutcut::Vector{Int8}
  cell_to_bgcell::Vector{Int32}
  point_to_coords::Vector{Point{D,T}}
  point_to_rcoords::Vector{Point{D,T}}
end

struct FacetSubTriangulation{D,T}
  facet_to_points::Table{Int}
  facet_to_normal::Vector{VectorValue{D,T}}
  facet_to_bgcell::Vector{Int32}
  point_to_coords::Vector{Point{D,T}}
  point_to_rcoords::Vector{Point{D,T}}
end

struct NewBulkMesh{D,T,M}
  background_mesh::M
  subtriangulation::SubTriangulation{D,T}
  facet_subtriangulations::Vector{FacetSubTriangulation{D,T}}
  background_cell_to_inoutcut::Vector{Int8}
end

struct BulkMesh{D,T,M}
  background_mesh::M
  cell_to_vertex_coordinates::Table{Point{D,T}}
  cell_to_d_to_dface_to_vertices::Table{Table{Int}}
  cell_to_d_to_dface_to_in_out_boundary::Vector{Table{Int}}
  cell_to_in_out_cut::Vector{Int}
  vertex_to_in_out_boundary::Vector{Int}
end

function NewBulkMesh(bg_mesh::M,sm::SurfaceMesh{D,T}) where {D,T,M}
  @check num_dims(bg_mesh) == D

  cell_mesh = CellMesh{D,T}()

  c_to_v = Table{Int}()
  c_to_io = Int8[]
  c_to_bgc = Int32[]
  c_coords = Point{D,T}[]


  f_to_v = Table{Int}()
  f_to_n = VectorValue{D,T}[]
  f_to_bgc = Int32[]
  f_coords = Point{D,T}[]

  bgc_to_ioc = fill(Int8(FACE_UNDEF),num_cells(bg_mesh))
  bgv_to_iob = fill(Int8(FACE_UNDEF),num_vertices(bg_mesh))

  c_to_sm_f = compute_cell_to_surface_mesh_faces(bg_mesh,sm)
  
  for k in 1:num_cells(bg_mesh)
    reset!(cell_mesh, get_cell(bg_mesh,k) )
    compute_cell_mesh!(cell_mesh,sm,c_to_sm_f,k)
    if is_surface_mesh_captured(cell_mesh)
     
      bgc_to_ioc[k] = FACE_CUT
      for lv in 1:num_vertices_per_cell(bg_mesh)
        v = get_vertex_id(bg_mesh,k,lv)
        if bgv_to_iob[v] == FACE_UNDEF
          bgv_to_iob[v] = get_vertex_in_out_boundary(cell_mesh,lv)
        end
      end
      
      offset = length(c_coords)
      append!( c_coords, get_vertex_coordinates(cell_mesh) )
      append!( c_to_v, get_cell_to_vertices(cell_mesh), offset )
      
      for cell in 1:num_cells(cell_mesh)
        push!( c_to_bgc, k )
        if is_cell_interior(cell_mesh,cell)
          push!(c_to_io, FACE_IN )
        elseif is_cell_exterior(cell_mesh,cell)
          push!(c_to_io, FACE_OUT )
        else
          throw(ErrorException(""))
        end
      end

      vertices_cache = Int[]
      for facet in 1:num_facets(cell_mesh)
        if is_facet_boundary(cell_mesh,facet)
          resize!(vertices_cache, 0 )
          facet_coordinates = get_facet_coordinates(cell_mesh,facet)
          
          push!( f_to_bgc, k ) 
          push!( f_to_n, normal(facet_coordinates) )
          for v in get_vertices(facet_coordinates)
            push!(f_coords,v)
            push!(vertices_cache, length(f_coords))
          end
          push!(f_to_v, vertices_cache )
        end
      end

    end

  end

  subtriangulation = SubTriangulation( c_to_v, c_to_io, c_to_bgc, c_coords, c_coords )
  facet_subtriangulation = FacetSubTriangulation( f_to_v, f_to_n, f_to_bgc, f_coords, f_coords )
  facet_subtriangulations = [ facet_subtriangulation ]
  
  ## Create ST's

  queue = Int[]
  cells_around = Int[]

  for cell in 1:num_cells(bg_mesh)
    if bgc_to_ioc[cell] == FACE_CUT
      head = 1
      resize!(queue,1)
      queue[head] = cell

      while head ≤ length(queue)
        current_cell = queue[head]
        head += 1
        for lvertex in 1:num_vertices_per_cell(bg_mesh)
          vertex = get_vertex_id(bg_mesh,current_cell,lvertex)
          if bgv_to_iob[vertex] == FACE_IN
            for cell_around in cells_around_vertex!(cells_around,bg_mesh,vertex)
              if bgc_to_ioc[cell_around] == FACE_UNDEF
                bgc_to_ioc[cell_around] = FACE_IN
                for _lvertex in 1:num_vertices_per_cell(bg_mesh)
                  _vertex = get_vertex_id(bg_mesh,cell_around,_lvertex)
                  if bgv_to_iob[_vertex] == FACE_UNDEF
                    bgv_to_iob[_vertex] = FACE_IN
                  end
                  push!(queue,cell_around)
                end
              end
            end
          end 
        end
      end
    end
  end

  for cell in 1:num_cells(bg_mesh)
    if bgc_to_ioc[cell] == FACE_UNDEF
      bgc_to_ioc[cell] = FACE_OUT
    end
  end

  NewBulkMesh{D,T,M}( bg_mesh, subtriangulation, facet_subtriangulations, bgc_to_ioc )

end



function BulkMesh(bg_mesh::M,sm::SurfaceMesh{D,T}) where {D,T,M}
  @check num_dims(bg_mesh) == D

  cell_mesh = CellMesh{D,T}()

  c_to_vc = Table{Point{D,T}}()
  c_to_d_to_df_to_v = Table{Table{Int}}()
  c_to_d_to_df_to_io = Table{Int}[]
  c_to_io = fill(FACE_UNDEF,num_cells(bg_mesh))
  v_to_io = fill(FACE_UNDEF,num_vertices(bg_mesh))

  c_to_sm_f = compute_cell_to_surface_mesh_faces(bg_mesh,sm)
  
  for k in 1:num_cells(bg_mesh)
    reset!(cell_mesh, get_cell(bg_mesh,k) )
    compute_cell_mesh!(cell_mesh,sm,c_to_sm_f,k)
    if is_surface_mesh_captured(cell_mesh)
     
      c_to_io[k] = FACE_CUT
      for lv in 1:num_vertices_per_cell(bg_mesh)
        v = get_vertex_id(bg_mesh,k,lv)
        if v_to_io[v] == FACE_UNDEF
          v_to_io[v] = get_vertex_in_out_boundary(cell_mesh,lv)
        end
      end
      
      vertices = copy(get_vertex_coordinates(cell_mesh))
      push!( c_to_vc, vertices )

      tables = [ Table{Int}() for d in 0:num_dims(cell_mesh) ]
      for d in 0:num_dims(cell_mesh)
        append!(tables[d+1],get_dface_to_vertices(cell_mesh,d))
      end
      push!( c_to_d_to_df_to_v, tables )

      table = deepcopy(cell_mesh.d_to_dface_to_in_out_boundary)
      push!( c_to_d_to_df_to_io, table )
      
    else
      push!( c_to_vc, [] )
      push!( c_to_d_to_df_to_v, [] )
      push!( c_to_d_to_df_to_io, Table{Int}() )
    end

  end

  queue = Int[]
  cells_around = Int[]

  for cell in 1:num_cells(bg_mesh)
    if c_to_io[cell] == FACE_CUT
      head = 1
      resize!(queue,1)
      queue[head] = cell

      while head ≤ length(queue)
        current_cell = queue[head]
        head += 1
        for lvertex in 1:num_vertices_per_cell(bg_mesh)
          vertex = get_vertex_id(bg_mesh,current_cell,lvertex)
          if v_to_io[vertex] == FACE_IN
            for cell_around in cells_around_vertex!(cells_around,bg_mesh,vertex)
              if c_to_io[cell_around] == FACE_UNDEF
                c_to_io[cell_around] = FACE_IN
                for _lvertex in 1:num_vertices_per_cell(bg_mesh)
                  _vertex = get_vertex_id(bg_mesh,cell_around,_lvertex)
                  if v_to_io[_vertex] == FACE_UNDEF
                    v_to_io[_vertex] = FACE_IN
                  end
                  push!(queue,cell_around)
                end
              end
            end
          end 
        end
      end
    end
  end

  for cell in 1:num_cells(bg_mesh)
    if c_to_io[cell] == FACE_UNDEF
      c_to_io[cell] = FACE_OUT
    end
  end

  BulkMesh( bg_mesh, c_to_vc, c_to_d_to_df_to_v, c_to_d_to_df_to_io, c_to_io, v_to_io )
end

num_dims(m::BulkMesh{D}) where D = D

get_background_mesh(m::BulkMesh) = m.background_mesh

get_background_mesh(m::NewBulkMesh) = m.background_mesh

num_cells(m::BulkMesh) = num_cells( get_background_mesh(m) )

num_cells(m::NewBulkMesh) = num_cells( get_background_mesh(m) )

num_vertices(m::BulkMesh) = num_vertices( get_background_mesh(m) )

num_vertices(m::NewBulkMesh) = num_vertices( get_background_mesh(m) )

num_vertices_per_cell(m::BulkMesh) = num_vertices_per_cell( get_background_mesh(m) )

num_vertices_per_cell(m::NewBulkMesh) = num_vertices_per_cell( get_background_mesh(m) )

num_subvertices(m::NewBulkMesh) = length( get_subvertex_coordinates(m) )

num_subcells(m::NewBulkMesh) = length( m.subtriangulation.cell_to_points ) 

get_subvertex_coordinates(m::NewBulkMesh) = m.subtriangulation.point_to_coords

get_subcell_to_in_out_cut(m::NewBulkMesh,cell::Integer) = m.subtriangulation.cell_to_inoutcut[cell]

function get_vertex_id(m::BulkMesh,cell::Integer,lvertex::Integer)
  bg_mesh = get_background_mesh(m)
  get_vertex_id(bg_mesh,cell,lvertex)
end

function get_vertex_id(m::NewBulkMesh,cell::Integer,lvertex::Integer)
  bg_mesh = get_background_mesh(m)
  get_vertex_id(bg_mesh,cell,lvertex)
end

get_cell_in_out(m::BulkMesh,cell::Integer) = m.cell_to_in_out_cut[cell]

is_cell_cut(m::BulkMesh,cell::Integer) = m.cell_to_in_out_cut[cell] == FACE_CUT

is_cell_cut(m::NewBulkMesh,cell::Integer) = m.background_cell_to_inoutcut[cell] == FACE_CUT

get_cell_to_in_out_cut(m::NewBulkMesh,cell::Integer) = m.background_cell_to_inoutcut[cell]

function get_vertex_coordinates(m::BulkMesh,i::Integer)
  bg_mesh = get_background_mesh(m)
  get_vertex_coordinates(bg_mesh,i)
end

function get_vertex_coordinates(m::NewBulkMesh,i::Integer)
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

function writevtk(m::NewBulkMesh{D,T},file_base_name) where {D,T}
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
  for (i,p) in enumerate( get_subvertex_coordinates(m) )
    for d in 1:D
      points[d,offset+i] = p[d]
    end
  end
  
  cells = MeshCell{Vector{Int64}}[]
  cell_to_io = Int8[]

  vtk_type = VTKCellType(d_to_vtk_hex_type_id[D])
  for cell in 1:num_cells(m)
    if !is_cell_cut(m,cell)
     vertices = zeros(Int,num_vertices_per_cell(m))
      for j in 1:num_vertices_per_cell(m)
        lv = vtk_hex_lvertices[j]
        vertices[j] = get_vertex_id(m,cell,lv)
      end
      push!( cells, MeshCell(vtk_type,vertices) )
      push!( cell_to_io, get_cell_to_in_out_cut(m,cell) )
    end
  end

  cell_to_vertices = m.subtriangulation.cell_to_points
  vtk_type = VTKCellType(d_to_vtk_tet_type_id[D])
  for cell in 1:num_subcells(m)
    vertices = [ offset+cell_to_vertices[cell,lvertex] for lvertex in 1:length(cell_to_vertices,cell) ]
    push!( cells, MeshCell(vtk_type,vertices) )
    push!( cell_to_io, get_subcell_to_in_out_cut(m,cell) )
  end
  
  vtkfile = vtk_grid(file_base_name,points,cells)
  vtkfile["IO",VTKCellData()] = cell_to_io
  vtk_save(vtkfile)


  writevtk(m.facet_subtriangulations[1],file_base_name*"_facets")
end


num_vertices(facets::FacetSubTriangulation) = length(facets.point_to_coords)

num_facets(facets::FacetSubTriangulation) = length(facets.facet_to_points)

get_vertex_coordinates(facets::FacetSubTriangulation) = facets.point_to_coords

get_normal(facets::FacetSubTriangulation,facet::Integer) = facets.facet_to_normal[facet]

function get_vertex(facets::FacetSubTriangulation,facet::Integer,lvertex::Integer)
  facets.facet_to_points[facet,lvertex]
end

function writevtk(facets::FacetSubTriangulation{D,T},file_base_name) where {D,T}
  d_to_vtk_tet_type_id = Dict(0=>1,1=>3,2=>5,3=>10)
  points = zeros(T,D,num_vertices(facets))
  normals = zeros(T,D,num_facets(facets))
  for (i,p) in enumerate( get_vertex_coordinates(facets) ), d in 1:D
    points[d,i] = p[d]
  end
  cells = MeshCell{Vector{Int64}}[]
  vtk_type = VTKCellType(d_to_vtk_tet_type_id[D-1])
  for facet in 1:num_facets(facets)
    vertices = [ get_vertex(facets,facet,lvertex) for lvertex in 1:vtk_type.nodes ]
    push!( cells, MeshCell(vtk_type,vertices) )
    n = get_normal(facets,facet)
    for d in 1:D
      normals[d,facet] = n[d]
    end
  end
  vtkfile = vtk_grid(file_base_name,points,cells)
  vtkfile["normals",VTKCellData()] = normals
  vtk_save(vtkfile)
end
