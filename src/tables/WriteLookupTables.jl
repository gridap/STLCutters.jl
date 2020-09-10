module WriteLookupTables

include("LookupTablesGenerator.jl")

LookupTablesGenerator.write_tables

dir = @__DIR__
filename = "LookupTables.jl"

filename = joinpath(dir,filename)

LookupTablesGenerator.write_tables(filename)

end # module
