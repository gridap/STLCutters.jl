module RunSaveTests

using Test
using STLCuttersAnalysis
using STLCuttersAnalysis: download_run_and_save
using DrWatson

filename = download_thing(80084,path=tmpdir())
data = run_and_save(filename,rerun=true,nmax=20,nmin=5)
@test data["surface_error"] < 1e-9
@test data["volume_error"] < 1e-9
@test !data["has_leak"]

filename = datadir("geometries","cube.stl")
data = run_and_save(filename,rerun=true,nmax=7,poisson=true)
@test data["error_l2"] < 1e-9
@test data["error_h1"] < 1e-9

end # module

