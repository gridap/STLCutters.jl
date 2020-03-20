
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

struct BulkMesh{D,T,M}
  background_mesh::M
  subtriangulation::SubTriangulation{D,T}
  facet_subtriangulations::Vector{FacetSubTriangulation{D,T}}
  background_cell_to_inoutcut::Vector{Int8}
end

function BulkMesh(bg_mesh::M,sm::SurfaceMesh{D,T}) where {D,T,M}
  @check num_dims(bg_mesh) == D

  cell_mesh = CellMesh{D,T}()

  c_to_v = Table{Int}()
  c_to_io = Int8[]
  c_to_bgc = Int32[]
  c_coords = Point{D,T}[]
  c_rcoords = Point{D,T}[]


  f_to_v = Table{Int}()
  f_to_n = VectorValue{D,T}[]
  f_to_bgc = Int32[]
  f_coords = Point{D,T}[]
  f_rcoords = Point{D,T}[]

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
      vertices = get_vertex_coordinates(cell_mesh) 
      ref_vertices = transformation( get_reference_cell(bg_mesh), get_cell(bg_mesh,k), vertices )
      append!( c_coords, vertices )
      append!( c_rcoords, ref_vertices )
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
          facet_normal = normal(facet_coordinates)
          facet_normal = facet_normal / norm(facet_normal)

          push!( f_to_bgc, k ) 

          push!( f_to_n, facet_normal )
          for v in get_vertices(facet_coordinates)
            push!(f_coords,v)
            push!(f_rcoords, transformation( get_reference_cell(bg_mesh), get_cell(bg_mesh,k), v ) )
            push!(vertices_cache, length(f_coords))
          end
          push!(f_to_v, vertices_cache )
        end
      end

    end

  end

  subtriangulation = SubTriangulation( c_to_v, c_to_io, c_to_bgc, c_coords, c_rcoords )
  facet_subtriangulation = FacetSubTriangulation( f_to_v, f_to_n, f_to_bgc, f_coords, f_rcoords )
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

  BulkMesh{D,T,M}( bg_mesh, subtriangulation, facet_subtriangulations, bgc_to_ioc )

end

get_background_mesh(m::BulkMesh) = m.background_mesh

num_cells(m::BulkMesh) = num_cells( get_background_mesh(m) )

num_vertices(m::BulkMesh) = num_vertices( get_background_mesh(m) )

num_cells(m::SubTriangulation) = length(m.cell_to_points)

get_vertex_coordinates(m::SubTriangulation) = m.point_to_coords

get_cell_to_vertices(m::SubTriangulation) = m.cell_to_points

get_cell_to_inout(m::SubTriangulation,cell::Integer) = m.cell_to_inoutcut[cell]

function get_cell_coordinates(m::SubTriangulation{D},cell::Integer) where D
  v =  get_vertex_coordinates(m)
  c_to_v = get_cell_to_vertices(m)
  @check length(c_to_v,1) == D+1 "get_cell_coordinates: only implemented for simplex mesh"
  if D == 2
    Triangle( v[c_to_v[cell,1]], v[c_to_v[cell,2]], v[c_to_v[cell,3]] )
  elseif D == 3
    Tetrahedron( v[c_to_v[cell,1]], v[c_to_v[cell,2]], v[c_to_v[cell,3]], v[c_to_v[cell,4]] )
  else
    throw(ErrorException("get_cell_coordinates(::SubTriangulation{$D}) not implemented"))
  end
end

num_vertices(facets::FacetSubTriangulation) = length(facets.point_to_coords)

num_facets(facets::FacetSubTriangulation) = length(facets.facet_to_points)

get_vertex_coordinates(facets::FacetSubTriangulation) = facets.point_to_coords

get_facet_to_vertices(facets::FacetSubTriangulation) = facets.facet_to_points

get_normal(facets::FacetSubTriangulation,facet::Integer) = facets.facet_to_normal[facet]

function get_vertex(facets::FacetSubTriangulation,facet::Integer,lvertex::Integer)
  facets.facet_to_points[facet,lvertex]
end

function get_facet_coordinates(facets::FacetSubTriangulation{D},cell::Integer) where D
  v =  get_vertex_coordinates(facets)
  f_to_v = get_facet_to_vertices(facets)
  @check length(f_to_v,1) == D  "get_facet_coordinates: only implemented for simplex mesh"
  if D == 2
    Segment( v[f_to_v[cell,1]], v[f_to_v[cell,2]] )
  elseif D == 3
    Triangle( v[f_to_v[cell,1]], v[f_to_v[cell,2]], v[f_to_v[cell,3]] )
  else
    throw(ErrorException("get_facet_coordinates(::SubTriangulation{$D}) not implemented"))
  end
end


function volume(m::SubTriangulation)
  volume = 0.0
  for cell in 1:num_cells(m)
    c_coords = get_cell_coordinates(m,cell)
    volume += measure(c_coords)
  end
  volume
end

function surface(facets::FacetSubTriangulation)
  surface = 0.0
  for facet in 1:num_facets(facets)
    f_coords = get_facet_coordinates(facets,facet)
    surface += measure(f_coords)
  end
  surface
end

surface(m::BulkMesh,i::Integer) = surface(m.facet_subtriangulations[i])

volume(m::BulkMesh) = volume(m.subtriangulation)

function writevtk(m::BulkMesh{D,T},file_base_name) where {D,T}
  for (i,facets) in enumerate( m.facet_subtriangulations )
    writevtk(facets,"$(file_base_name)_facets_$i")
  end
  
  ( points, cells, cell_to_io ) = get_vtk_cells( m.background_mesh, m.background_cell_to_inoutcut )
  ( subpoints, subcells, subcell_to_io ) = get_vtk_cells( m.subtriangulation, size(points,2) )
  points = hcat( points, subpoints )
  append!( cells, subcells )
  append!( cell_to_io, subcell_to_io )
  
  vtkfile = vtk_grid(file_base_name,points,cells)
  vtkfile["IO",VTKCellData()] = cell_to_io
  vtk_save(vtkfile)
end

function get_vtk_cells(m::SubTriangulation{D,T},offset::Integer) where {D,T}
  d_to_vtk_tet_type_id = Dict(0=>1,1=>3,2=>5,3=>10)
  points = zeros(Float32,D,num_cells(m))
  for (i,p) in enumerate( get_vertex_coordinates(m) )
    for d in 1:D
      points[d,i] = p[d]
    end
  end
  cells = MeshCell{Vector{Int64}}[]
  cell_to_io = Int8[]
  cell_to_vertices = m.cell_to_points
  vtk_type = VTKCellType(d_to_vtk_tet_type_id[D])
  for cell in 1:num_cells(m)
    vertices = [ offset+cell_to_vertices[cell,lvertex] for lvertex in 1:length(cell_to_vertices,cell) ]
    push!( cells, MeshCell(vtk_type,vertices) )
    push!( cell_to_io, get_cell_to_inout(m,cell) )
  end
  ( points, cells, cell_to_io )
end

function get_vtk_cells(m::CartesianMesh{D,T},cell_to_ioc::Vector) where {D,T}
  d_to_vtk_hex_type_id = Dict(0=>1,1=>3,2=>9,3=>12)
  vtk_hex_lvertices = [1,2,4,3,5,6,8,7]
  points = zeros(Float32,D,num_vertices(m))
  for i in 1:num_vertices(m)
    p = get_vertex_coordinates(m,i)
    for d in 1:D
      points[d,i] = p[d]
    end
  end
  cells = MeshCell{Vector{Int64}}[]
  cell_to_io = Int8[]
  vtk_type = VTKCellType(d_to_vtk_hex_type_id[D])
  for cell in 1:num_cells(m)
    if cell_to_ioc[cell] != FACE_CUT
     vertices = zeros(Int,num_vertices_per_cell(m))
      for j in 1:num_vertices_per_cell(m)
        lv = vtk_hex_lvertices[j]
        vertices[j] = get_vertex_id(m,cell,lv)
      end
      push!( cells, MeshCell(vtk_type,vertices) )
      push!( cell_to_io, cell_to_ioc[cell] )
    end
  end
  ( points, cells, cell_to_io )
end

function writevtk(facets::FacetSubTriangulation{D,T},file_base_name) where {D,T}
  d_to_vtk_tet_type_id = Dict(0=>1,1=>3,2=>5,3=>10)
  points = zeros(T,D,num_vertices(facets))
  normals = zeros(T,3,num_facets(facets))
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
