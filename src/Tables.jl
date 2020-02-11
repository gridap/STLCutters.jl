
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

Base.length(a::TableOfVectors) = length(a._vectors)

pushlist!(a::TableOfVectors{T},list::Vector{T}) where{T} = push!(a._vectors,list)

push_to_list!(a::TableOfVectors{T},i::Int,value::T) where{T} = push!(a._vectors[i],value)

set_to_list!(a::TableOfVectors{T},i::Int,j::Int,value::T) where {T} = a._vectors[i][j]=value

@inline Base.:(==)(a::TableOfVectors{T},b::TableOfVectors{T}) where T = ( a._vectors == b._vectors )

struct TableOfLists{T} <: AbstractTable
  _l::Vector{T}
  _p::Vector{Int}
end

function TableOfLists(x::Vector{Tuple{Int64,T}},n::Int) where T
  l = zeros(T,length(x))
  p = zeros(Int,n+1)
  for (i,v) in x
    p[i+1] += 1
  end
  p[1]=1
  for i in 1:n
    p[i+1] = p[i+1]+p[i]
  end
  for (i,v) in x
    l[p[i]] = v
    p[i] += 1
  end
  for i in 1:n-1
    p[n-i+1]=p[n-i]
  end
  p[1]=1
  TableOfLists(l,p)
end

table_cache(a::TableOfLists) = nothing

Base.length(a::TableOfLists) = length(a._p)-1

function STLCutter.getlist!(::Nothing,a::TableOfLists,i::Integer)
  a._l[a._p[i]:a._p[i+1]-1]
end
