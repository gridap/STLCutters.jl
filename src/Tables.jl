
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

  function Table(data::Vector{Vector{T}}) where T
    _data, _ptrs = compress_data(data)
    _masks = fill(true,length(data))
    new{T}(_data,_ptrs,_masks)
  end

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

function Base.push!(a::Table,b::Vector)
  push!(a.data,b)
  push!(a.masks,true)
  push!(a.ptrs,length(data)+1)
end

function remove!(a::Table,i::Integer)
  @check isactive(a,i)
  a.masks[i] = false
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


abstract type AbstractTable end

function table_cache(a::AbstractTable)
  @abstractmethod
end

function getlist!(cache,a::AbstractTable,i::Integer)
  @abstractmethod
end

function Base.length(a::AbstractTable)
  @abstractmethod
end

function getlist(a::AbstractTable,i::Integer)
  cache = table_cache(a)
  getlist!(cache,a,i)
end

function Base.show(io::IO,a::AbstractTable)
  print(io,"[")
  cache = table_cache(a)
  for i in 1:length(a)
    if i >1
      print(io,", ")
    end
    sub_list = getlist!(cache,a,i)
    print(io,sub_list)
  end
  print(io,"]")
end

struct TableOfVectors{T} <: AbstractTable
  _vectors::Vector{Vector{T}}
end

table_cache(a::TableOfVectors) = nothing

function getlist!(::Nothing,a::TableOfVectors,i::Integer)
  a._vectors[i]
end

Base.getindex(a::TableOfVectors,i::Int,j::Int) = a._vectors[i][j]

Base.setindex!(a::TableOfVectors{T},i::Int,j::Int,v::T) where T = a._vectors[i][j]=v

Base.length(a::TableOfVectors) = length(a._vectors)

Base.length(a::TableOfVectors,i::Int) = length(a._vectors[i])

@inline Base.:(==)(a::TableOfVectors{T},b::TableOfVectors{T}) where T = ( a._vectors == b._vectors )

