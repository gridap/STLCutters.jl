module BulkMeshesTests

using STLCutter

using STLCutter: expand, FACE_UNDEF, reset!, add_surface_mesh_face!,compute_in_out!,is_surface_mesh_captured, get_vertex_coordinates, BulkMesh, FACE_CUT, num_vertices_per_cell, get_vertex_id, get_vertex_in_out_boundary, FACE_IN, cells_around_vertex!, FACE_OUT, @check

using Test


stl = STL(joinpath(@__DIR__,"data/cube.stl"))
#stl = STL(joinpath(@__DIR__,"data/Bunny-LowPoly.stl"))

sm = SurfaceMesh(stl)

box = BoundingBox(sm)

box = expand(box,0.5)

bg_mesh = CartesianMesh(box,9)
cell_mesh = CellMesh(box)

c_to_sm_f = compute_cell_to_surface_mesh_faces(bg_mesh,sm)

D = num_dims(bg_mesh)
T = Float64

c_to_vc = Table{Point{D,T}}()
c_to_d_to_df_to_v = Table{Table{Int}}()
c_to_d_to_df_to_io = Table{Int}[]
c_to_io = fill(FACE_UNDEF,num_cells(bg_mesh))
v_to_io = fill(FACE_UNDEF,num_vertices(bg_mesh))

for k in 1:num_cells(bg_mesh)
  cell_coordinates = get_cell(bg_mesh,k)
  reset!(cell_mesh,cell_coordinates)
  for i in 1:length(c_to_sm_f,k)
    sm_face = c_to_sm_f[k,i]
    add_surface_mesh_face!(cell_mesh,sm,sm_face)
  end
  compact!(cell_mesh)
  compute_in_out!(cell_mesh,sm)

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

bm = BulkMesh( bg_mesh, c_to_vc, c_to_d_to_df_to_v, c_to_d_to_df_to_io, c_to_io, v_to_io )

writevtk(bm,"bulk")

end # module
