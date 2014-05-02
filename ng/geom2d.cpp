#include "stdafx.h" 
#include "mystdlib.h"

#include "myadt.hpp"
#include "gprim.hpp"

namespace netgen
{

ostream & operator<<(ostream  & s, const Point2d & p)
{
  return s << "(" << p.px << ", " << p.py << ")";
}

ostream & operator<<(ostream  & s, const Vec2d & v)
{
  return s << "(" << v.vx << ", " << v.vy << ")";
}

double Fastatan2 (double x, double y)
{
  if (y > 0)
    {
      if (x > 0)
	return y / (x+y);
      else
	return 1 - x / (y-x);
    }
  else if (y < 0)
    {
      if (x < 0)
	return 2 + y / (x+y);
      else
	return 3 - x / (y-x);
    }
  else 
    {
      if (x >= 0)
	return 0;
      else
	return 2;
    }
}

double Angle (const Vec2d & v)
{
  if (v.X() == 0 && v.Y() == 0)
    return 0;
    
  double ang = atan2 (v.Y(), v.X());
  if (ang < 0) ang+= 2 * M_PI;
  return ang;
}

double FastAngle (const Vec2d & v)
{
  return Fastatan2 (v.X(), v.Y());
}

double Angle (const Vec2d & v1, const Vec2d & v2)
{
  double ang = Angle(v2) - Angle(v1);
  if (ang < 0) ang += 2 * M_PI;
  return ang;
}

double FastAngle (const Vec2d & v1, const Vec2d & v2)
{
  double ang = FastAngle(v2) - FastAngle(v1);
  if (ang < 0) ang += 4;
  return ang;
}

double Dist2(const Line2d & g, const Line2d & h )
  {
  double   dd = 0.0, d1,d2,d3,d4;
  Point2d  cp = CrossPoint(g,h);
  
  if ( Parallel(g,h) || !IsOnLine(g,cp) || !IsOnLine(h,cp) )
    {
      d1 = Dist2(g.P1(),h.P1());
      d2 = Dist2(g.P1(),h.P2());
      d3 = Dist2(g.P2(),h.P1());
      d4 = Dist2(g.P2(),h.P2());
      if (d1<d2)  d2 = d1;
      if (d3<d4)  d4 = d3;
      dd = ( d2 < d4 ) ? d2 : d4;
    }
  return dd;
}

Point2d CrossPoint (const Line2d & l1, const Line2d & l2)
  {
  double den = Cross (l1.Delta(), l2.Delta());
  double num = Cross ( (l2.P1() - l1.P1()), l2.Delta());

  if (den == 0)
    return l1.P1();
  else
    return l1.P1() + (num/den) * l1.Delta();
}

int CrossPointBarycentric (const Line2d & l1, const Line2d & l2,
			   double & lam1, double & lam2)
{
  // p = l1.1 + lam1 (l1.2-l1.1) = l2.1 + lam2 (l2.2-l2.1)
  double a11 = l1.p2.X() - l1.p1.X();
  double a21 = l1.p2.Y() - l1.p1.Y();
  double a12 = -(l2.p2.X() - l2.p1.X());
  double a22 = -(l2.p2.Y() - l2.p1.Y());

  double b1 = l2.p1.X() - l1.p1.X();
  double b2 = l2.p1.Y() - l1.p1.Y();
  
  double det = a11*a22 - a12 * a21;
  if (det == 0)
    return 1;
    
  lam1 = (a22 * b1 - a12 * b2) / det;
  lam2 = (a11 * b2 - a21 * b1) / det;
  return 0;
}

int Parallel (const Line2d & l1, const Line2d & l2, double peps)
  {
  double p = fabs (Cross (l1.Delta(), l2.Delta()));
  return p <= peps * l1.Length() * l2.Length();
}

int IsOnLine (const Line2d & l, const Point2d & p, double heps)
  {
  double c1 = (p - l.P1()) * l.Delta();
  double c2 = (p - l.P2()) * l.Delta();
  double d = fabs (Cross ( (p - l.P1()), l.Delta()));
  double len2 = l.Length2();

  return c1 >= -heps * len2 && c2 <= heps * len2 && d <= heps * len2;
}

Polygon2d :: Polygon2d ()
{
  ;
}

void Polygon2d :: AddPoint (const Point2d & p)
{ 
  points.Append(p); 
}

double Polygon2d :: HArea () const
{
  int i;
  double ar = 0;
  for (i = 1; i <= points.Size(); i++)
    {
      const Point2d & p1 = points.Get(i);
      const Point2d & p2 = points.Get(i%points.Size()+1);
      ar += 
	(p2.X()-p1.X()) * p1.Y() -
	(p2.Y()-p1.Y()) * p1.X();
    }
  return ar/2;
}

int Polygon2d :: IsOn (const Point2d & p) const
{
  int i;
  for (i = 1; i <= points.Size(); i++)
    {
      const Point2d & p1 = points.Get(i);
      const Point2d & p2 = points.Get(i%points.Size()+1);
      if (IsOnLine (Line2d(p1, p2), p)) return 1;
    }
  return 0;
}

int Polygon2d :: IsIn (const Point2d & p) const
{
  int i;
  double sum = 0, ang;
  for (i = 1; i <= points.Size(); i++)
    {
      const Point2d & p1 = points.Get(i);
      const Point2d & p2 = points.Get(i%points.Size()+1);
      ang = Angle ( (p1 - p), (p2 - p) );
      if (ang > M_PI) ang -= 2 * M_PI;
      sum += ang;
    }
  return fabs(sum) > M_PI;
}

Point2d Polygon2d :: Center () const
  {
    return Point2d (0, 0);
  }
}