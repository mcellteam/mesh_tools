BASEDIR = ..
INC = -I${BASEDIR}/include

CPP = g++-3.3
#CPP = g++
#CPPFLAGS = ${INC} -m64 -O3 -march=opteron -felide-constructors -Wall -Wpointer-arith -Wwrite-strings -Wold-style-cast
#CPPFLAGS = ${INC} -O3 -march=pentium4  -felide-constructors -Wall -Wpointer-arith -Wwrite-strings -Wold-style-cast
CPPFLAGS = ${INC} -O3 -felide-constructors -Wall -Wpointer-arith -Wwrite-strings -Wold-style-cast
#CPPFLAGS = ${INC} -O3 -march=pentium4  -felide-constructors -Wall -Wpointer-arith -Wwrite-strings -Wold-style-cast
#CPPFLAGS = ${INC} -g -march=pentium4  -Wall -Wpointer-arith -Wwrite-strings -Wold-style-cast
#CFLAGS = -O3 -march=pentium4 -mcpu=pentium4
CFLAGS = -O3
#CFLAGS = -g -march=pentium4 -mcpu=pentium4
LDFLAGS = -static
HLIBS = libHh.a ../lib/recipes/librecipes.a ../lib/linpack/liblinpack.a
LIBS = -lm

.SUFFIXES: .C .c .o
.C.o:
	$(CPP) $(CPPFLAGS) -c $<

PROGS =	Filtermesh Filtera3d Recon Meshfit Subdivfit Polyfit
# order is important: Hh.o Timer.o last so their pools empty early.
LIBHh =	\
	Stack.o Queue.o Pqueue.o Map.o Set.o \
	Combination.o HashStruct.o UnionFind.o Graph.o GraphOp.o \
	Pool.o Stat.o Random.o \
	Geometry.o Matrix4.o Homogeneous.o GeomOp.o \
	Polygon.o Bbox.o \
	HashFloat.o HashPoint.o Spatial.o PolygonSpatial.o \
	GMesh.o Mesh.o MeshOp.o Facedistance.o SubMesh.o \
	A3dStream.o FrameIO.o BA3dStream.o Buffer.o \
	Mk3d.o Mklib.o \
	Options.o FileIO.o \
	LLS.o Principal.o Recipes.o \
	Kdtree.o Contour.o Hh.o Timer.o

all:	${PROGS}
libHh.a: ${LIBHh}
	ar r $@ ${LIBHh}
	ranlib $@
clean:
	rm -f ${PROGS} *.o *.a
OBFiltera3d =	Filtera3d.o ${HLIBS}
Filtera3d:	${OBFiltera3d}
	$(CPP) $(CPPFLAGS) $(LDFLAGS) ${OBFiltera3d} -o $@ ${LIBS}
OBFiltermesh =	Filtermesh.o ${HLIBS}
Filtermesh:	${OBFiltermesh}
	$(CPP) $(CPPFLAGS) $(LDFLAGS) ${OBFiltermesh} -o $@ ${LIBS}
OBRecon =	Recon.o ${HLIBS}
Recon:	${OBRecon}
	$(CPP) $(CPPFLAGS) $(LDFLAGS) ${OBRecon} -o $@ ${LIBS}
OBMeshfit =	Meshfit.o MFFit.o MFProject.o MFLocalFit.o ${HLIBS}
Meshfit:	${OBMeshfit}
	$(CPP) $(CPPFLAGS) $(LDFLAGS) ${OBMeshfit} -o $@ ${LIBS}
OBPolyfit =	Polyfit.o ${HLIBS}
Polyfit:	${OBPolyfit}
	$(CPP) $(CPPFLAGS) $(LDFLAGS) ${OBPolyfit} -o $@ ${LIBS}
OBSubdivfit =	Subdivfit.o ${HLIBS}
Subdivfit:	${OBSubdivfit}
	$(CPP) $(CPPFLAGS) $(LDFLAGS) ${OBSubdivfit} -o $@ ${LIBS}

# fast matrix inverse
Matrix4.o:	Matrix4.C
	$(CPP) $(CPPFLAGS) -c $<
# fast g3d hlr
Kdtree.o:	Kdtree.C
	$(CPP) $(CPPFLAGS) -c $<
# for ProjectPTri() **
Facedistance.o:	Facedistance.C
	$(CPP) $(CPPFLAGS) -c $<
# for ProjectPNeighb
MeshOp.o:	MeshOp.C
	$(CPP) $(CPPFLAGS) -c $<
# for global fit
LLS.o:	LLS.C
	$(CPP) $(CPPFLAGS) -c $<
# for LocalFit()
MFLocalFit.o:	MFLocalFit.C
	$(CPP) $(CPPFLAGS) -c $<
# for ProjectPointFaces()
MFProject.o:	MFProject.C
	$(CPP) $(CPPFLAGS) -c $<
# Pqueue very heavy in Meshfit
Pqueue.o:	Pqueue.C
	$(CPP) $(CPPFLAGS) -c $<
# for ObjectSpatial::enter(), ObjectSpatial::addcell()
Spatial.o:	Spatial.C
	$(CPP) $(CPPFLAGS) -c $<

# DO NOT DELETE THIS LINE -- g++dep uses it.
# DO NOT PUT ANYTHING AFTER THIS LINE, IT WILL GO AWAY.

A3dStream.o : A3dStream.C \
  Hh.h \
  A3dStream.h \
  Geometry.h \
  Array.h \
  Pool.h
BA3dStream.o : BA3dStream.C \
  Hh.h \
  BA3dStream.h \
  A3dStream.h \
  Geometry.h \
  Array.h \
  Pool.h \
  Buffer.h
Bbox.o : Bbox.C \
  Hh.h \
  Bbox.h \
  Geometry.h \
  Array.h \
  Pool.h
Buffer.o : Buffer.C \
  Hh.h \
  Buffer.h
Combination.o : Combination.C \
  Hh.h \
  Combination.h \
  Map.h \
  Pool.h \
  Stack.h
Contour.o : Contour.C \
  Hh.h \
  Contour.h \
  Geometry.h \
  Array.h \
  Pool.h \
  GMesh.h \
  Mesh.h \
  Map.h \
  Set.h \
  Stack.h \
  Pqueue.h \
  Polygon.h \
  Queue.h \
  MeshOp.h \
  Stat.h
Facedistance.o : Facedistance.C \
  Hh.h \
  Facedistance.h \
  Geometry.h \
  Array.h \
  Pool.h
FileIO.o : FileIO.C \
  Hh.h \
  FileIO.h
Filtera3d.o : Filtera3d.C \
  Hh.h \
  Options.h \
  Queue.h \
  Pool.h \
  FileIO.h \
  A3dStream.h \
  Geometry.h \
  Array.h \
  FrameIO.h \
  Stat.h \
  Bbox.h \
  Random.h \
  Polygon.h \
  Kdtree.h \
  HashPoint.h \
  HashFloat.h \
  Map.h \
  HashStruct.h \
  Graph.h \
  Stack.h \
  Spatial.h \
  Set.h \
  Pqueue.h
Filterframe.o : Filterframe.C \
  Hh.h \
  Options.h \
  Queue.h \
  Pool.h \
  FrameIO.h \
  Geometry.h \
  Array.h
Filtermesh.o : Filtermesh.C \
  Hh.h \
  Options.h \
  Queue.h \
  Pool.h \
  GMesh.h \
  Mesh.h \
  Map.h \
  Set.h \
  Stack.h \
  Pqueue.h \
  Geometry.h \
  Array.h \
  MeshOp.h \
  Stat.h \
  A3dStream.h \
  FileIO.h \
  HashPoint.h \
  HashFloat.h \
  HashStruct.h \
  Homogeneous.h \
  Matrix4.h \
  GeomOp.h \
  Polygon.h \
  Random.h \
  Timer.h \
  Principal.h
FrameIO.o : FrameIO.C \
  Hh.h \
  FrameIO.h \
  Geometry.h \
  Array.h \
  Pool.h \
  Buffer.h
GMesh.o : GMesh.C \
  Hh.h \
  GMesh.h \
  Mesh.h \
  Map.h \
  Pool.h \
  Set.h \
  Stack.h \
  Pqueue.h \
  Geometry.h \
  Array.h \
  Polygon.h \
  A3dStream.h
GeomOp.o : GeomOp.C \
  Hh.h \
  GeomOp.h \
  Geometry.h \
  Array.h \
  Pool.h
Geometry.o : Geometry.C \
  Hh.h \
  Geometry.h \
  Array.h \
  Pool.h \
  Homogeneous.h \
  Matrix4.h
Graph.o : Graph.C \
  Hh.h \
  Graph.h \
  Map.h \
  Pool.h \
  Stack.h
GraphOp.o : GraphOp.C \
  Hh.h \
  Geometry.h \
  Array.h \
  Pool.h \
  GraphOp.h \
  Graph.h \
  Map.h \
  Stack.h \
  Pqueue.h \
  Set.h \
  UnionFind.h \
  Spatial.h \
  Stat.h \
  Queue.h
HX.o : HX.C \
  Hh.h \
  HX.h \
  Options.h \
  Queue.h \
  Pool.h \
  HX.xbm
HashFloat.o : HashFloat.C \
  Hh.h \
  HashFloat.h \
  Map.h \
  Pool.h
HashPoint.o : HashPoint.C \
  Hh.h \
  HashPoint.h \
  HashFloat.h \
  Map.h \
  Pool.h \
  HashStruct.h \
  Geometry.h \
  Array.h
HashStruct.o : HashStruct.C \
  Hh.h \
  HashStruct.h \
  Map.h \
  Pool.h
Hh.o : Hh.C \
  Hh.h \
  Map.h \
  Pool.h
Hlr.o : Hlr.C \
  Hh.h \
  Hlr.h \
  Geometry.h \
  Array.h \
  Pool.h \
  Polygon.h \
  Kdtree.h \
  Bbox.h
Homogeneous.o : Homogeneous.C \
  Hh.h \
  Homogeneous.h \
  Geometry.h \
  Array.h \
  Pool.h \
  Matrix4.h
Kdtree.o : Kdtree.C \
  Hh.h \
  Kdtree.h \
  Array.h \
  Stack.h \
  Pool.h \
  Stat.h
LLS.o : LLS.C \
  Hh.h \
  LLS.h \
  Array.h \
  Stack.h \
  Pool.h \
  Recipes.h \
  ../include/nr.h \
  ../include/nrutil.h \
  Linpack.h \
  ../include/linpack.h \
  Timer.h
MFFit.o : MFFit.C \
  Hh.h \
  Meshfit.h \
  GMesh.h \
  Mesh.h \
  Map.h \
  Pool.h \
  Set.h \
  Stack.h \
  Pqueue.h \
  Geometry.h \
  Array.h \
  LLS.h
MFLocalFit.o : MFLocalFit.C \
  Hh.h \
  Meshfit.h \
  GMesh.h \
  Mesh.h \
  Map.h \
  Pool.h \
  Set.h \
  Stack.h \
  Pqueue.h \
  Geometry.h \
  Array.h \
  Facedistance.h \
  GeomOp.h \
  Stat.h
MFProject.o : MFProject.C \
  Hh.h \
  Meshfit.h \
  GMesh.h \
  Mesh.h \
  Map.h \
  Pool.h \
  Set.h \
  Stack.h \
  Pqueue.h \
  Geometry.h \
  Array.h \
  PolygonSpatial.h \
  Spatial.h \
  Polygon.h \
  Facedistance.h \
  MeshOp.h \
  Stat.h \
  Queue.h
Map.o : Map.C \
  Hh.h \
  Map.h \
  Pool.h \
  Set.h \
  Random.h
Matrix4.o : Matrix4.C \
  Hh.h \
  Matrix4.h
Mesh.o : Mesh.C \
  Hh.h \
  Mesh.h \
  Map.h \
  Pool.h \
  Set.h \
  Stack.h \
  Pqueue.h \
  Geometry.h \
  Array.h \
  Random.h
MeshOp.o : MeshOp.C \
  Hh.h \
  MeshOp.h \
  GMesh.h \
  Mesh.h \
  Map.h \
  Pool.h \
  Set.h \
  Stack.h \
  Pqueue.h \
  Geometry.h \
  Array.h \
  Stat.h \
  Queue.h \
  GeomOp.h \
  HashStruct.h \
  Polygon.h \
  Homogeneous.h \
  Matrix4.h \
  Facedistance.h
Meshfit.o : Meshfit.C \
  Hh.h \
  Meshfit.h \
  GMesh.h \
  Mesh.h \
  Map.h \
  Pool.h \
  Set.h \
  Stack.h \
  Pqueue.h \
  Geometry.h \
  Array.h \
  Options.h \
  Queue.h \
  MeshOp.h \
  Stat.h \
  A3dStream.h \
  FileIO.h \
  GeomOp.h \
  Bbox.h \
  Facedistance.h \
  FrameIO.h \
  HashStruct.h \
  Random.h \
  Recipes.h \
  ../include/nr.h \
  ../include/nrutil.h \
  Timer.h
Mk3d.o : Mk3d.C \
  Hh.h \
  Mk3d.h \
  A3dStream.h \
  Geometry.h \
  Array.h \
  Pool.h \
  Stack.h
Mklib.o : Mklib.C \
  Hh.h \
  Mklib.h \
  Mk3d.h \
  A3dStream.h \
  Geometry.h \
  Array.h \
  Pool.h \
  Stack.h
Options.o : Options.C \
  Hh.h \
  Options.h \
  Queue.h \
  Pool.h
Polyfit.o : Polyfit.C \
  Hh.h \
  Options.h \
  Queue.h \
  Pool.h \
  A3dStream.h \
  Geometry.h \
  Array.h \
  FileIO.h \
  Bbox.h \
  FrameIO.h \
  LLS.h \
  Stack.h \
  Principal.h \
  Polygon.h \
  Set.h \
  Map.h \
  Random.h \
  Timer.h
Polygon.o : Polygon.C \
  Hh.h \
  Polygon.h \
  Geometry.h \
  Array.h \
  Pool.h \
  A3dStream.h \
  Bbox.h
PolygonSpatial.o : PolygonSpatial.C \
  Hh.h \
  PolygonSpatial.h \
  Spatial.h \
  Geometry.h \
  Array.h \
  Pool.h \
  Map.h \
  Stack.h \
  Set.h \
  Pqueue.h \
  Polygon.h \
  Facedistance.h \
  Bbox.h
Pool.o : Pool.C \
  Hh.h \
  Pool.h
Postscript.o : Postscript.C \
  Hh.h \
  Postscript.h \
  Geometry.h \
  Array.h \
  Pool.h
Pqueue.o : Pqueue.C \
  Hh.h \
  Pqueue.h \
  Map.h \
  Pool.h
Principal.o : Principal.C \
  Hh.h \
  Principal.h \
  Geometry.h \
  Array.h \
  Pool.h \
  Homogeneous.h \
  Matrix4.h \
  Recipes.h \
  ../include/nr.h \
  ../include/nrutil.h
Queue.o : Queue.C \
  Hh.h \
  Queue.h \
  Pool.h
Random.o : Random.C \
  Hh.h \
  Random.h
Recipes.o : Recipes.C \
  Hh.h \
  Recipes.h \
  ../include/nr.h \
  ../include/nrutil.h
Recon.o : Recon.C \
  Hh.h \
  Options.h \
  Queue.h \
  Pool.h \
  A3dStream.h \
  Geometry.h \
  Array.h \
  Homogeneous.h \
  Matrix4.h \
  Mk3d.h \
  Stack.h \
  Mklib.h \
  Stat.h \
  Polygon.h \
  Spatial.h \
  Map.h \
  Set.h \
  Pqueue.h \
  Principal.h \
  Graph.h \
  GraphOp.h \
  Contour.h \
  GMesh.h \
  Mesh.h \
  MeshOp.h \
  Timer.h \
  Bbox.h \
  FrameIO.h
Repository.o : Repository.C \
  Hh.h
Set.o : Set.C \
  Hh.h \
  Set.h \
  Map.h \
  Pool.h
Spatial.o : Spatial.C \
  Hh.h \
  Spatial.h \
  Geometry.h \
  Array.h \
  Pool.h \
  Map.h \
  Stack.h \
  Set.h \
  Pqueue.h \
  Bbox.h \
  Queue.h \
  Stat.h
Stack.o : Stack.C \
  Hh.h \
  Stack.h \
  Pool.h
Stat.o : Stat.C \
  Hh.h \
  Stat.h \
  Queue.h \
  Pool.h
SubMesh.o : SubMesh.C \
  Hh.h \
  SubMesh.h \
  GMesh.h \
  Mesh.h \
  Map.h \
  Pool.h \
  Set.h \
  Stack.h \
  Pqueue.h \
  Geometry.h \
  Array.h \
  Combination.h \
  Homogeneous.h \
  Matrix4.h \
  MeshOp.h \
  Stat.h \
  Queue.h \
  HashStruct.h
Subdivfit.o : Subdivfit.C \
  Hh.h \
  Options.h \
  Queue.h \
  Pool.h \
  A3dStream.h \
  Geometry.h \
  Array.h \
  FileIO.h \
  PolygonSpatial.h \
  Spatial.h \
  Map.h \
  Stack.h \
  Set.h \
  Pqueue.h \
  Polygon.h \
  SubMesh.h \
  GMesh.h \
  Mesh.h \
  Combination.h \
  Homogeneous.h \
  Matrix4.h \
  Facedistance.h \
  MeshOp.h \
  Stat.h \
  GeomOp.h \
  LLS.h \
  Timer.h \
  HashStruct.h \
  Recipes.h \
  ../include/nr.h \
  ../include/nrutil.h \
  Random.h \
  Bbox.h \
  FrameIO.h
Template.o : Template.C \
  Hh.h \
  Map.h \
  Pool.h
Timer.o : Timer.C \
  Hh.h \
  Timer.h \
  Map.h \
  Pool.h \
  Stat.h \
  Queue.h
UnionFind.o : UnionFind.C \
  Hh.h \
  UnionFind.h \
  Map.h \
  Pool.h

# IF YOU PUT ANYTHING HERE IT WILL GO AWAY
