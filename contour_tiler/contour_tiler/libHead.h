#ifndef LIB_HEAD_H
#define LIB_HEAD_H

class HexOutput{
	public:
		int* hexes;
		int* faces;
		float* verts;

		int numHexes;
		int numFaces;
		int numVerts;
};

class TileConfig {
};

HexOutput* getHexas(int no_of_tris, int* tri, int no_of_pts, float* pnts);

#endif