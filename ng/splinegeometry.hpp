#pragma once

#include "spline.hpp"

namespace netgen
{


  template < int D >
  class SplineGeometry 
  {
    // protected:
  public:  
    Array < GeomPoint<D> > geompoints;
    Array < SplineSeg<D>* > splines;

     ~SplineGeometry();

     int Load (const Array<double> & raw_data, const int startpos = 0);

    const Array<SplineSeg<D>*> & GetSplines () const
    { return splines; }

    int GetNSplines (void) const { return splines.Size(); }
    string GetSplineType (const int i) const { return splines[i]->GetType(); }
    SplineSeg<D> & GetSpline (const int i) {return *splines[i];}
    const SplineSeg<D> & GetSpline (const int i) const {return *splines[i];}

     void GetBoundingBox (Box<D> & box) const;
    Box<D> GetBoundingBox () const 
    { Box<D> box; GetBoundingBox (box); return box; }

    int GetNP () const { return geompoints.Size(); }
    const GeomPoint<D> & GetPoint(int i) const { return geompoints[i]; }

    void AppendSegment(SplineSeg<D> * spline)
    {
      splines.Append (spline);
    }
  };

}