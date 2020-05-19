
struct Table{T}
  data::Vector{T}
  ptrs::Vector{Int32}
  masks::Vector{Bool}

  function Table(data::Vector,rows::Vector{Int},n::Int)
    T = eltype(data)
    _data, _ptrs = compress_data(data,rows,n)
    _masks = fill(true,n)
    new{T}(_data,_ptrs,_masks)
  end

  function Table(data::Vector,ptrs::Vector)
    T = eltype(data)
    masks = fill(true,length(ptrs)-1)
    new{T}(data,ptrs,masks)
  end

  function Table{T}(data::Vector{Vector{T}}) where T
    _data, _ptrs = compress_data(data)
    _masks = fill(true,length(_ptrs)-1)
    new{T}(_data,_ptrs,_masks)
  end

  Table(data::Vector{Vector{T}}) where T = Table{T}(data)
  
  function Table(data::Matrix)
    T = eltype(data)
    _data, _ptrs = compress_data(data)
    _masks = fill(true,length(_ptrs)-1)
    new{T}(_data,_ptrs,_masks)
  end

end

function Table{T}() where T
  _data = T[]
  _ptrs = Int32[1]
  Table(_data,_ptrs)
end

Table(::Type{T}) where T = Table{T}() 

function Table{T}(::UndefInitializer,n::Integer,m::Integer) where T
  ptrs = fill(Int32(m),n+1)
  length_to_ptrs!(ptrs)
  data = Vector{T}(undef,ptrs[n+1]-1)
  Table(data,ptrs)
end

Base.length(a::Table) = length(a.ptrs)-1

@inline Base.length(a::Table,i::Integer) = a.ptrs[i+1] - a.ptrs[i]

Base.eltype(::Type{Table{T}}) where T = T

Base.eltype(::Table{T}) where T = T

@inline function Base.getindex(a::Table,i::Integer,j::Integer)
  @check isactive(a,i)
  p = a.ptrs[i]-1
  a.data[p+j]
end

@inline function Base.setindex!(a::Table,v,i::Integer,j::Integer)
  @check isactive(a,i)
  p = a.ptrs[i]-1
  a.data[p+j] = v
end

function Base.show(io::IO,a::Table)
  T = eltype(a)
  print(io,"$(typeof(a))( Vector{$T}[")
  for i in 1:length(a)
    if isactive(a,i)
      i == 1 || print(io,",")
      print(io," [")
      for j in 1:length(a,i)
        j == 1 || print(io,",")
        print(io," $(a[i,j])")
      end
      print(io," ]")
    end
  end
  print(io," ] )")
end

function Base.show(io::IO,::MIME"text/plain",a::Table)
  T = typeof(a)
  n = length(a)
  println(io,"$n-row $T:")
  for i in 1:n
    print(io,"[")
    if isactive(a,i)
      for j in 1:length(a,i)
        print(io,a[i,j])
        if j < length(a,i)
          print(io,",")
        end
      end
    else
      print(io,"-")
    end
    println(io,"]")
  end
end

function Base.:(==)(a::Table,b::Table)
  (a.data == b.data) && (a.ptrs == b.ptrs) && (a.masks == b.masks)
end
 
Base.maximum(a::Table) = maximum(a.data)

function Base.resize!(a::Table,n::Int)
  @check n <= length(a)
  resize!(a.ptrs,n+1)
  resize!(a.masks,n)
  resize!(a.data,a.ptrs[n+1]-1)
  a
end

function Base.fill!(a::Table,val)
  fill!(a.data,val)
  a
end

function Base.resize!(a::Table,len::Vector)
  n = length(a)
  resize!(a.ptrs,1)
  append!(a.ptrs,len)
  length_to_ptrs!(a.ptrs)
  resize!(a.data,a.ptrs[end]-1)
  resize!(a.masks,a.ptrs[end]-1)
  fill!(a.masks,true)
  a
end

function Base.push!(a::Table,b::Vector)
  append!(a.data,b)
  push!(a.masks,true)
  push!(a.ptrs,length(a.data)+1)
  a
end

function Base.append!(a::Table,b::Matrix)
  n = size(b,2)
  m = size(b,1)
  for i in 1:n
    for j in 1:m
      push!(a.data,b[j,i])
    end
    push!(a.ptrs,length(a.data)+1)
    push!(a.masks,true)
  end
  a
end

function Base.append!(a::Table,b::Table)
  for i in 1:length(b)
    if isactive(b,i)
      for j in 1:length(b,i)
        push!(a.data,b[i,j])
      end
      push!(a.ptrs,length(a.data)+1)
      push!(a.masks,true)
    end
  end
  a
end

function Base.append!(a::Table,b::Table,offset)
  for i in 1:length(b)
    if isactive(b,i)
      for j in 1:length(b,i)
        push!(a.data,b[i,j]+offset)
      end
      push!(a.ptrs,length(a.data)+1)
      push!(a.masks,true)
    end
  end
  a
end

function Base.append!(a::Table,b::Vector{<:Vector})
  for i in 1:length(b)
    push!(a,b[i])
  end
  a
end

function remove!(a::Table,i::Integer)
  @check isactive(a,i)
  a.masks[i] = false
  a
end

@inline function isactive(a::Table,i::Integer)
  a.masks[i]
end

function compact!(a::Table)
  k = 0
  for i in 1:length(a)
    if isactive(a,i)
      for j in 1:length(a,i)
        k +=1
        a.data[k] = a[i,j]
      end
    end
  end
  resize!(a.data,k)
  k = 0
  for i in 1:(length(a.ptrs)-1)
    if isactive(a,i)
      k +=1
      n = a.ptrs[i+1] - a.ptrs[i]
      a.ptrs[k+1] = a.ptrs[k] + n
    end
  end
  resize!(a.ptrs,k+1)
  resize!(a.masks,k)
  fill!(a.masks,true)
end

# Helpers

function compress_data(data::Vector,rows::Vector{Int},n::Int)
  @check length(data) == length(rows)
  T = eltype(data)
  _ptrs = zeros(Int32,n+1)
  for r in rows
    _ptrs[r+1] +=  1
  end
  length_to_ptrs!(_ptrs)
  ndata = _ptrs[end]-1
  _data = zeros(T,ndata)
  for k in 1:length(data)
    r = rows[k]
    p = _ptrs[r]
    _data[p] = data[k]
    _ptrs[r] += 1
  end
  rewind_ptrs!(_ptrs)
  (_data, _ptrs)
end

function rewind_ptrs!(ptrs::AbstractVector{<:Integer})
  @inbounds for i in (length(ptrs)-1):-1:1
    ptrs[i+1] = ptrs[i]
  end
  ptrs[1] = 1
end

function length_to_ptrs!(ptrs::AbstractArray{<:Integer})
  ptrs[1] = 1
  @inbounds for i in 1:(length(ptrs)-1)
    ptrs[i+1] += ptrs[i]
  end
end

function compress_data(data::Vector{Vector{T}}) where T
  n = length(data)
  _ptrs = zeros(Int32,n+1)
  for i in 1:n
    _ptrs[i+1] = length(data[i])
  end
  length_to_ptrs!(_ptrs)
  ndata = _ptrs[end]-1
  _data = zeros(T,ndata)
  c = 0
  for i in 1:n, j in 1:length(data[i])
    _data[c+=1] = data[i][j]
  end
  (_data,_ptrs)
end

function compress_data(data::Matrix)
  n = size(data,2)
  m = size(data,1)
  _ptrs = fill(Int32(m),n+1)
  length_to_ptrs!(_ptrs)
  _data = [ data[i] for i in 1:length(data) ]
  (_data,_ptrs)
end
