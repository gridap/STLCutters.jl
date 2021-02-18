
analysisdir(args...) = projectdir("analysis",args...)

testdir(args...) = projectdir("test",args...)

DrWatson.datadir(args...) = analysisdir("data",args...)

DrWatson.scriptsdir(args...) = analysisdir("scripts",args...)

DrWatson.plotsdir(args...) = analysisdir("plots",args...)

tmpdir(args...) = analysisdir("tmp",args...)

