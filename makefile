
all:
	cd contour_tiler/contour_tiler; make
	cd contour_tiler/poly_converter; make
	cd dx2mesh; make
	cd dx2mesh_2; make
	cd irit2mesh; make
	cd mesh2dx; make
	cd mesh2irit; make
	cd mesh2mcell; make
	cd mesh2off; make
	cd mesh2rib; make
	cd mesh2ribwireframe; make
	cd mesh2smesh; make
	cd mesh2smf; make
	cd mesh2stl; make
	cd mesh2vtk; make
	cd meshclip; make
	cd meshfilter; make
	cd meshmerge; make
	cd meshoffset; make
	cd meshscale; make
	cd meshtranslate; make
	cd netgen2mesh; make
	cd netgen2smesh; make
	cd Recon_tmb_1; make
	cd simplify_surface; make
	cd stl2mesh; make
	cd stlb2stla; make
	cd synu2mesh; make
	cd vizvtk; make
	cd vrml2mesh; make
	cd vtk2mesh; make
clean:
	cd contour_tiler/contour_tiler; make clean
	cd contour_tiler/poly_converter; make clean
	cd dx2mesh; make clean
	cd dx2mesh_2; make clean
	cd irit2mesh; make clean
	cd mesh2dx; make clean
	cd mesh2irit; make clean
	cd mesh2mcell; make clean
	cd mesh2off; make clean
	cd mesh2rib; make clean
	cd mesh2ribwireframe; make clean
	cd mesh2smesh; make clean
	cd mesh2smf; make clean
	cd mesh2stl; make clean
	cd mesh2vtk; make clean
	cd meshclip; make clean
	cd meshfilter; make clean
	cd meshmerge; make clean
	cd meshoffset; make clean
	cd meshscale; make clean
	cd meshtranslate; make clean
	cd netgen2mesh; make clean
	cd netgen2smesh; make clean
	cd Recon_tmb_1; make clean
	cd simplify_surface; make clean
	cd stl2mesh; make clean
	cd stlb2stla; make clean
	cd synu2mesh; make clean
	cd vizvtk; make clean
	cd vrml2mesh; make clean
	cd vtk2mesh; make clean
