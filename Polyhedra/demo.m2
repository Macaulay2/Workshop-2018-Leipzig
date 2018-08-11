restart
path
path = prepend("~/GitRepos/PolyhedraWork/M2/Macaulay2/packages", path)
loadPackage "Polyhedra"

viewHelp Polyhedra

-- Some basic functionality
P = convexHull transpose matrix {{0,0},{2,0},{2,3},{-1,4}}
vertices P
facets P
latticePoints P
isEmpty P
isCompact P
F = normalFan P
rays F
isSimplicial F
isSimplicial P
isSmooth F

-- 'cellDecompose' has become 'regularSubdivision'
V = transpose matrix {{1,1},{0,0},{1,0},{0,1}}
w = matrix{{1,1,0,0}}
L1 = regularSubdivision(convexHull V, w)
vertices L1#0
L2 = regularSubdivision(V, w)
-- Why are these different?
latticePoints convexHull V
-- These have a different order than the points in V
w = matrix{{0,1,0,1}}
regularSubdivision(convexHull V, w)
regularSubdivision(V, w)
-- This method was previously called cellDecompose
cellDecompose(V, w)
cellDecompose(convexHull V, w)

-- 'triangulate' has become 'regularTriangulation' and 'barycentricTriangulation' (adds vertices)
S = convexHull V
triangulate S
barycentricTriangulation S
regularTriangulation S -- calls TOPCOM

-- added 'minimalNonFaces' method
F = normalFan hypercube 3
rays F
minimalNonFaces F
stanleyReisnerRing F

-- added parameters for hypercubes
vertices hypercube 2
vertices hypercube(2,0,1) -- (dimension, lower-bound, upper-bound)
vertices hypercube(2,-1,0)

-- added many tests (now > 400)
-- added documentation

