#ifndef VORONOICURVE_H_
#define VORONOICURVE_H_

// standard includes
#include <iostream>
#include <fstream>
#include <cassert>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// includes for defining the Voronoi diagram adaptor
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include "Voronoi_diagram_2.h"

#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>


//#include <GL/glut.h>
#include <stdlib.h>
#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif



namespace Voronoicrv
{

typedef CGAL::Exact_predicates_inexact_constructions_kernel                  K;
typedef CGAL::Delaunay_triangulation_2<K>                                    DT;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>                 AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
typedef CGAL::Voronoi_diagram_2<DT,AT,AP>                                    VD;

// typedef for the result type of the point location
typedef AT::Site_2                    Site_2;
typedef AT::Point_2                   Point_2;

typedef VD::Locate_result             Locate_result;
typedef VD::Vertex_handle             Vertex_handle;
typedef VD::Face_handle               Face_handle;
typedef VD::Halfedge_handle           Halfedge_handle;
typedef VD::Ccb_halfedge_circulator   Ccb_halfedge_circulator;

typedef VD::Face_iterator Face_iterator;
typedef VD::Edge_iterator Edge_iterator;
typedef VD::Vertex_iterator Vertex_iterator;

using namespace std;

class VoronoiCurve
{
private:
	VD vd;
	vector<pair<int,int>> _boundary;
        vector<pair<double, double> > points;
	void reconstruct();
	void collectBoundaryIndices(int pci, int cind, Point_2 peelpcon[][2], Point_2 con[][2]);
	int getIndex(double, double);
        int iterator1(VD::Vertex_iterator vit);
        VD::Vertex_iterator find_it(Point_2 p);

public:
	VoronoiCurve(vector<pair<double, double> > pointVec);
	vector<pair<int,int>> *getBoundary();
};
}

#endif /* VORONOICURVE_H_ */
