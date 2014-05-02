#include "stdafx.h" 
#include "mystdlib.h"

#include "csg.hpp"
#include "meshing.hpp"

namespace netgen
{

Meshing2Surfaces :: Meshing2Surfaces (const Surface & asurf,
				      const MeshingParameters & mp,
				      const Box<3> & abb)
  : Meshing2(mp, abb), surface(asurf)
{
  ;
}

void Meshing2Surfaces :: DefineTransformation (const Point3d & p1, const Point3d & p2,
					       const PointGeomInfo * geominfo1,
					       const PointGeomInfo * geominfo2)
{
  ((Surface&)surface).DefineTangentialPlane (p1, p2);
}

void Meshing2Surfaces :: TransformToPlain (const Point3d & locpoint, 
					   const MultiPointGeomInfo & geominfo,
					   Point2d & planepoint, 
					   double h, int & zone)
{
  Point<2> hp;
  surface.ToPlane (locpoint, hp, h, zone);
  planepoint.X() = hp(0);
  planepoint.Y() = hp(1);
}

int Meshing2Surfaces :: TransformFromPlain (Point2d & planepoint,
					    Point3d & locpoint, 
					    PointGeomInfo & gi,
					    double h)
{
  Point<3> hp;
  Point<2> hp2 (planepoint.X(), planepoint.Y());
  surface.FromPlane (hp2, hp, h);
  locpoint = hp;
  gi.trignum = 1;
  return 0;
}

double Meshing2Surfaces :: CalcLocalH (const Point3d & p, double gh) const
{
  return surface.LocH (p, 3, 1, gh);
}
}