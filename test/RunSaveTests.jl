module RunSaveTests

using Test
using STLCutters
using STLCutters.Tests
using STLCutters.Tests: download_run_and_save


filename = joinpath(@__DIR__,"data/47076.stl")
data = run_and_save(filename,rerun=true,nmax=50)
@test data["surface_error"] < 1e-9
@test data["volume_error"] < 1e-9
@test !data["has_leak"]

filename = joinpath(@__DIR__,"data/47076_sf.obj")
data = run_and_save(filename,rerun=true,nmax=50)
@test data["surface_error"] < 1e-9
@test data["volume_error"] < 1e-9
@test !data["has_leak"]

filename = download_thing(293137)
data = run_and_save(filename,rerun=true,nmax=50)
@test data["surface_error"] < 1e-9
@test data["volume_error"] < 1e-9
@test !data["has_leak"]


filename = download_thing(80084)
data = run_and_save(filename,rerun=true,nmax=20,nmin=5)
@test data["surface_error"] < 1e-9
@test data["volume_error"] < 1e-9
@test !data["has_leak"]

filename = download_thing(65395)
data = run_and_save(filename,rerun=true,nmax=20,nmin=5)
@test data["surface_error"] < 1e-9
@test data["volume_error"] < 1e-9
@test !data["has_leak"]

filename = download_thing(77343)
data = run_and_save(filename,rerun=true,nmax=20,nmin=5)
@test data["surface_error"] < 1e-9
@test data["volume_error"] < 1e-9
@test !data["has_leak"]


filename = download_thing(95436)
data = run_and_save(filename,rerun=true,nmax=20,nmin=5)
@test data["surface_error"] < 1e-9
@test data["volume_error"] < 1e-9
@test !data["has_leak"]


filename = download_thing(93703)
data = run_and_save(filename,rerun=true,nmax=20,nmin=5)
@test data["surface_error"] < 1e-9
@test data["volume_error"] < 1e-9
@test !data["has_leak"]

filename = download_thing(243015)
data = run_and_save(filename,rerun=true,nmax=100,nmin=1)
@test data["surface_error"] < 1e-9
@test data["volume_error"] < 1e-9
@test !data["has_leak"]


end # module
