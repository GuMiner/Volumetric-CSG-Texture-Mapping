#include "stdafx.h" 
#include "mystdlib.h"
#include "myadt.hpp"

#include "linalg.hpp"
#include "gprim.hpp"

namespace netgen
{
void
LocalCoordinates (const Vec3d & e1, const Vec3d & e2,
		  const Vec3d & v, double & lam1, double & lam2)
{
  double m11 = e1 * e1;
  double m12 = e1 * e2;
  double m22 = e2 * e2;
  double rs1 = v * e1;
  double rs2 = v * e2;
  
  double det = m11 * m22 - m12 * m12;
  lam1 = (rs1 * m22 - rs2 * m12)/det;
  lam2 = (m11 * rs2 - m12 * rs1)/det;
}

double MinDistLP2 (const Point2d & lp1, const Point2d & lp2, const Point2d & p)
{
  Vec2d v(lp1, lp2);
  Vec2d vlp(lp1, p);

  // dist(lam) = \| vlp \|^2 - 2 lam (v1p, v) + lam^2 \| v \|^2

  // lam = (v * vlp) / (v * v);
  // if (lam < 0) lam = 0;
  // if (lam > 1) lam = 1;

  double num = v*vlp;
  double den = v*v;

  if (num <= 0) 
    return Dist2 (lp1, p);

  if (num >= den) 
    return Dist2 (lp2, p);
  
  if (den > 0)
    {
      return vlp.Length2() - num * num /den;
    }
  else
    return vlp.Length2();
}

double MinDistLP2 (const Point3d & lp1, const Point3d & lp2, const Point3d & p)
{
  Vec3d v(lp1, lp2);
  Vec3d vlp(lp1, p);

  // dist(lam) = \| vlp \|^2 - 2 lam (v1p, v) + lam^2 \| v \|^2

  // lam = (v * vlp) / (v * v);
  // if (lam < 0) lam = 0;
  // if (lam > 1) lam = 1;

  double num = v*vlp;
  double den = v*v;

  if (num <= 0) 
    return Dist2 (lp1, p);

  if (num >= den) 
    return Dist2 (lp2, p);
  
  if (den > 0)
    {
      return vlp.Length2() - num * num /den;
    }
  else
    return vlp.Length2();
}

double MinDistTP2 (const Point3d & tp1, const Point3d & tp2, 
		   const Point3d & tp3, const Point3d & p)
{
  double lam1, lam2;
  double res;

  LocalCoordinates (Vec3d (tp1, tp2), Vec3d (tp1, tp3),
		    Vec3d (tp1, p), lam1, lam2);
  int in1 = lam1 >= 0;
  int in2 = lam2 >= 0;
  int in3 = lam1+lam2 <= 1;
  
  if (in1 && in2 && in3)
    {
      Point3d pp = tp1 + lam1 * Vec3d(tp1, tp2) + lam2 *  Vec3d (tp1, tp3);
      res = Dist2 (p, pp);
    }
  else
    {
      res = Dist2 (tp1, p);
      if (!in1)
	{
	  double hv = MinDistLP2 (tp1, tp3, p);
	  if (hv < res) res = hv; 
	}
      if (!in2)
	{
	  double hv = MinDistLP2 (tp1, tp2, p);
	  if (hv < res) res = hv; 
	}
      if (!in3)
	{
	  double hv = MinDistLP2 (tp2, tp3, p);
	  if (hv < res) res = hv; 
	}
    }

  return res;

  Vec3d pp1(tp1, p);
  Vec3d v1(tp1, tp2), v2(tp1, tp3);

  double c = pp1.Length2();
  double cx = -2 * (pp1 * v1);
  double cy = -2 * (pp1 * v2);
  double cxx = v1.Length2();
  double cxy = 2 * (v1 * v2);
  double cyy = v2.Length2();

  QuadraticPolynomial2V pol (-c, -cx, -cy, -cxx, -cxy, -cyy);
  double res2 =  - pol.MaxUnitTriangle ();

  if (fabs (res - res2) > 1e-8)
    cout << "res and res2 differ: " << res << " != " << res2 << endl;
  return res2;
}

}