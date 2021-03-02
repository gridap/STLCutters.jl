

function intersection(h::Face{D,D},Π::AbstractVector{<:Plane},inout,tol) where D
  T = intersection(h,first(Π),inout,tol)
  if length(Π) == 1
    Tr = T
  else
    Tr = empty_mesh(T)
    Πn = view(Π,2:length(Π))
    for i in 1:length(T[1])
      k = Cell(T[1][i],[h.vertex_coordinates;T[2]],T[4])
      io = T[3][i]
      Tk = intersection(k,Πn,io,tol)
      append_mesh!(Tr,Tk)
    end
  end
  Tr
end

function intersection(k::Face{D,D},Π::Plane,inout,tol) where D
  marching_cell(k,Π,inout,tol)
end

function marching_cell(k::Face{D,D},Π::Plane,inout::Integer,tol) where D
  @notimplementedif !is_n_cube(k)
  case = compute_case(k,Π)
  Tnew = compute_new_cells(k,case)
  if ( case == 1 || case == (1<<(1<<D))+1 ) && inout ≠ UNSET
    k_to_io = fill( inout, length(Tnew) )
  else
    k_to_io = collect(get_cell_to_inout_from_case(D,case))
  end
  Xnew,Xnew_to_old = compute_new_vertices!(Tnew,k,Π,case,tol)
  Tnew,Xnew,k_to_io,get_polytope(k)
end

function compute_case(k::Face{D,D},Π::Plane) where D
  case = 0
  for i in 1:num_vertices(k)
    if signed_distance(k[i],Π) < 0
      case |= (1<<(i-1))
    end
  end
  case + 1
end

function compute_new_cells(k::Cell{D,D},case) where D
  compute_new_cells(k.nodes,k.vertex_coordinates,k.polytope,case)
end

function compute_new_cells(k::Cell{D,D},case) where D
  compute_new_cells(k.nodes,k.vertex_coordinates,k.polytope,case)
end

function compute_new_vertices!(T,k::Face{D,D},plane,case,tol) where D
  compute_new_vertices!(T,k.nodes,k.vertex_coordinates,k.polytope,plane,case,atol=tol)
end

function empty_mesh(a::Tuple)
  T,X,T_io,p = a 
  empty(T),X,empty(T_io),p
end

function append_mesh!(a::T,b::T) where T<:Tuple
  @assert a[end] == b[end]
  append!(a[1],b[1])
  append!(a[2],b[2])
  append!(a[3],b[3])
  a
end
