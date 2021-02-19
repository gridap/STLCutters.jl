using Mustache
using DrWatson

@quickactivate "STLCutters"

using STLCutters
using STLCutters.Tests


stl_list = [
  "cube",
  "wine_glass",
  "37881", # Gear
  "47076", # Azteca pyramid
  "96457", # NightHawk
  "550964", # Eiffel Tower
  "293137", # Bunny low-poly
  "441708", # Standford Bunny
  "35269", # Octocat
  "65904", # Heart 
  "551021", # Art de Triomph
  "37266", # Extruded Earth
  "252119"] # Angel

δ = 0.2
nmin = 1
nmax = 14 .* 2 .^ (0:5) # Matching bounding box with cells 1+2δ 

i = 2
filename = joinpath(testdir("data"),"$(stl_list[i]).stl")
rotations_and_displacements(
  filename,
  nmax=28,
  nmin=1,
  displacements=exp10.(-17:-10),
  angles=exp10.(-17:-10))

