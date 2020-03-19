
macro check(test)
  quote
    @boundscheck @assert $(esc(test)) $(string(test))
  end
end

macro check(test,msg)
  quote
    @boundscheck @assert $(esc(test)) $msg
  end
end

macro abstractmethod()
  quote
    error("This function belongs to an interface definition and cannot be used.")
  end
end

"""
    tfill(v, ::Val{D}) where D
Returns a tuple of length `D` that contains `D` times the object `v`.
In contrast to `tuple(fill(v,D)...)` which returns the same result, this function is type-stable.
"""
function tfill(v, ::Val{D}) where D
  t = tfill(v, Val{D-1}())
  (v,t...)
end

tfill(v,::Val{0}) = ()
tfill(v,::Val{1}) = (v,)
tfill(v,::Val{2}) = (v,v)
tfill(v,::Val{3}) = (v,v,v)

