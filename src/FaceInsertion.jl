# Symbolic n-cube refinement

function vertex_refinement(k,Vk) # octree_refinement
  # Implement the cell mesh Tk defined as the octree refinement rule with interior vertex v, subcells are defined by its Morton index and the vertex coordinates per level
  # We could consider a look-up table that provides two different
  # situations: 1) do nothing if vertex on vertex (even though
  # this part should have been processed before); 2) split in two
  # parts if vertex on edge; 3) splint in four parts in vertex on # face; 4) split in 8 parts if vertex in cell... easily to write
  # for a general dim
  return Tk
end

function edge_refinement(k,Ek)
  # look up table for all cases possible:
  # edge cut cells on boundary by construction
  # end-points m1 and m2
  # look up table for all cases
  # we could consider a partition of the cell into four cells
  # and in some situations some of these cells will be zero and
  # thus neglected. All this can be done automatically by using all
  # possible combination of m1 m2 owners... e.g.
  # m1 in v1 , m2 in v2, etc (putting always m1 and m2 in the middle of the
  # vef), e.g., if m1 in f1 and m2 in f2, consider m1 in mid-point of f1 and
  # m2 in mid-point of m2... This way, given a polytope, we can automatically
  # generate a partition into degenerated cubes
  # we could do the same for vertices... pretty simple, thus, in the tree there
  # will be non-active/empty/undef cells when they are logically empty
  #
  # In this process, we can not only provide a partition of the cell but 
  # also all facets too, using same ideas recursively for all dims
  return Tk
end

function face_refinement(k,Fk)
  # idem, the symbolic look up table generator is absolutely general!
  # the plane slicing cube problem has 6 different cases in the look up table
  # we must be able to recover each of those.
  # Finally, we can also define the singular cases in which intersection vertices
  # are cube vertices
  # https://pdfs.semanticscholar.org/be26/7e00110d096745fd9cff5d59791df6322126.pdf
  # E.g., consider that the node is now in the mid-point of one of the 
  # edges containing it and eliminate the resulting zero volume meshes (all
  # symbolic here)
  # For these six cases, we can easily create degenerated hex meshes (even by
  # hand)
end

# ...

function d_face_refinement(d,k,Fk)
  # some 4D meshes comments on the internet
  # https://stackoverflow.com/questions/44939879/how-should-i-handle-morphing-4d-objects-in-opengl/44970550#44970550
  # it seems some people have done this before at least in 4D
end







