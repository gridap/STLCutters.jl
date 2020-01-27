
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

table_cache(a::TableOfVectors) = nothing

function getlist!(::Nothing,a::TableOfVectors,i::Integer)
  a._vectors[i]
end

Base.length(a::TableOfVectors) = length(a._vectors)
