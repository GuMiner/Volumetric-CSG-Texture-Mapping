#pragma once
namespace netgen
{

  ///
  class Meshing2Surfaces : public Meshing2
  {
    ///
    const Surface & surface;
  
  public:
    ///
    //  Meshing2Surfaces (const Surface & asurf);
    ///
    Meshing2Surfaces (const Surface & asurf, const MeshingParameters & mp, 
		      const Box<3> & aboundingbox);

  protected:
    ///
    virtual void DefineTransformation (const Point3d & p1, const Point3d & p2,
				       const PointGeomInfo * geominfo1,
				       const PointGeomInfo * geominfo2);
    ///
    virtual void TransformToPlain (const Point3d & locpoint, 
				   const MultiPointGeomInfo & geominfo,
				   Point2d & plainpoint, 
				   double h, int & zone);
    ///
    virtual int TransformFromPlain (Point2d & plainpoint,
				    Point3d & locpoint, 
				    PointGeomInfo & gi,
				    double h);
    ///
    virtual double CalcLocalH (const Point3d & p, double gh) const;
  };
}