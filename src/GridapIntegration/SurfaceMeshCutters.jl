
struct SurfaceMeshCutter <: GridapEmbedded.Cutter end

function cut(cutter::SurfaceMeshCutter,background::DiscreteModel,geom::STLGeometry)
  data = #...
  EmbeddedDiscretization(background, data..., geom)
end

function cut(background::DiscreteModel,geom::STLGeometry)
  cutter = SurfaceMeshCutter()
  cut(cutter,background,geom)
end


