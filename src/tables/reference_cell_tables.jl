#
# Reference TET & HEX cells 
#  Please do not modify it by hand:
#    code automatically generated by table_generator.jl
#

const D_to_reference_coordinates_for_tet = Matrix{Float64}[
  [0.0 1.0] , 
 
  [0.0 1.0 0.0; 0.0 0.0 1.0] , 
 
  [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0] , 
 
  [0.0 1.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0 1.0] , 
]

const D_to_d_to_dface_to_vertices_for_tet = Vector{Vector{Vector{Int64}}}[
  [] , 
 
  [[[1, 2],[1, 3],[2, 3],],] , 
 
  [[[1, 2],[1, 3],[2, 3],[1, 4],[2, 4],[3, 4],],[[1, 2, 3],[1, 2, 4],[1, 3, 4],[2, 3, 4],],] , 
 
  [[[1, 2],[1, 3],[2, 3],[1, 4],[2, 4],[3, 4],[1, 5],[2, 5],[3, 5],[4, 5],],[[1, 2, 3],[1, 2, 4],[1, 3, 4],[2, 3, 4],[1, 2, 5],[1, 3, 5],[2, 3, 5],[1, 4, 5],[2, 4, 5],[3, 4, 5],],[[1, 2, 3, 4],[1, 2, 3, 5],[1, 2, 4, 5],[1, 3, 4, 5],[2, 3, 4, 5],],] , 
]

const D_to_reference_coordinates_for_hex = Matrix{Float64}[
  [-1.0 1.0] , 
 
  [-1.0 1.0 -1.0 1.0; -1.0 -1.0 1.0 1.0] , 
 
  [-1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0; -1.0 -1.0 1.0 1.0 -1.0 -1.0 1.0 1.0; -1.0 -1.0 -1.0 -1.0 1.0 1.0 1.0 1.0] , 
 
  [-1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0; -1.0 -1.0 1.0 1.0 -1.0 -1.0 1.0 1.0 -1.0 -1.0 1.0 1.0 -1.0 -1.0 1.0 1.0; -1.0 -1.0 -1.0 -1.0 1.0 1.0 1.0 1.0 -1.0 -1.0 -1.0 -1.0 1.0 1.0 1.0 1.0; -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0] , 
]

const D_to_d_to_dface_to_vertices_for_hex = Vector{Vector{Vector{Int64}}}[
  [] , 
 
  [[[1, 2],[3, 4],[1, 3],[2, 4],],] , 
 
  [[[1, 2],[3, 4],[1, 3],[2, 4],[5, 6],[7, 8],[5, 7],[6, 8],[1, 5],[2, 6],[3, 7],[4, 8],],[[1, 2, 3, 4],[5, 6, 7, 8],[1, 2, 5, 6],[3, 4, 7, 8],[1, 3, 5, 7],[2, 4, 6, 8],],] , 
 
  [[[1, 2],[3, 4],[1, 3],[2, 4],[5, 6],[7, 8],[5, 7],[6, 8],[1, 5],[2, 6],[3, 7],[4, 8],[9, 10],[11, 12],[9, 11],[10, 12],[13, 14],[15, 16],[13, 15],[14, 16],[9, 13],[10, 14],[11, 15],[12, 16],[1, 9],[2, 10],[3, 11],[4, 12],[5, 13],[6, 14],[7, 15],[8, 16],],[[1, 2, 3, 4],[5, 6, 7, 8],[1, 2, 5, 6],[3, 4, 7, 8],[1, 3, 5, 7],[2, 4, 6, 8],[9, 10, 11, 12],[13, 14, 15, 16],[9, 10, 13, 14],[11, 12, 15, 16],[9, 11, 13, 15],[10, 12, 14, 16],[1, 2, 9, 10],[3, 4, 11, 12],[1, 3, 9, 11],[2, 4, 10, 12],[5, 6, 13, 14],[7, 8, 15, 16],[5, 7, 13, 15],[6, 8, 14, 16],[1, 5, 9, 13],[2, 6, 10, 14],[3, 7, 11, 15],[4, 8, 12, 16],],[[1, 2, 3, 4, 5, 6, 7, 8],[9, 10, 11, 12, 13, 14, 15, 16],[1, 2, 3, 4, 9, 10, 11, 12],[5, 6, 7, 8, 13, 14, 15, 16],[1, 2, 5, 6, 9, 10, 13, 14],[3, 4, 7, 8, 11, 12, 15, 16],[1, 3, 5, 7, 9, 11, 13, 15],[2, 4, 6, 8, 10, 12, 14, 16],],] , 
]

