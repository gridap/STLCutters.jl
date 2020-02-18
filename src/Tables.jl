
macro abstractmethod()
  quote
    error("This function belongs to an interface definition and cannot be used.")
  end
end

abstract type AbstractTable end

function table_cache(a::AbstractTable)
  @abstractmethod
end

function getlist!(cache,a::AbstractTable,i::Integer)
  @abstractmethod
end

function pushlist!(a::AbstractTable,list::Vector)
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

function TableOfVectors(::Type{T},ni::Integer,nj::Integer) where T
  data = [ zeros(T,nj) for i in 1:ni ]
  TableOfVectors{T}( data )
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

struct CompressedTable{T} <: AbstractTable
  data::Vector{T}
  ptrs::Vector{Int}
end

function CompressedTable(indices::Vector{Int},values::Vector{T},n::Int) where T
  if length(indices) != length(values) 
    throw(ArgumentError("CompressedTable: indices and values must have the same length"))
  end
  l = zeros(T,length(values))
  p = zeros(Int,n+1)
  for i in indices 
    p[i+1] += 1
  end
  p[1]=1
  for i in 1:n
    p[i+1] = p[i+1]+p[i]
  end
  for (j,id) in enumerate(indices) 
    l[p[id]] = values[j]
    p[id] += 1
  end
  for i in 1:n-1
    p[n-i+1]=p[n-i]
  end
  p[1]=1
  CompressedTable(l,p)
end


Base.length(a::CompressedTable) = length(a.ptrs)-1

Base.length(a::CompressedTable,i::Int) = a.ptrs[i+1] - a.ptrs[i]

Base.getindex(t::CompressedTable,i::Int,j::Int) = t.data[t.ptrs[i]+j-1]
