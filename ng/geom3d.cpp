#include "stdafx.h" 
#include "mystdlib.h"

#include <algorithm>
#include "myadt.hpp"
#include "gprim.hpp"

namespace netgen
{
ostream & operator<<(ostream  & s, const Point3d & p)
  {
  return s << "(" << p.x[0] << ", " << p.x[1] << ", " << p.x[2] << ")";
  }

ostream & operator<<(ostream  & s, const Vec3d & v)
  {
  return s << "(" << v.x[0] << ", " << v.x[1] << ", " << v.x[2] << ")";
  }

double Angle (const Vec3d & v1, const Vec3d & v2)
{
  double co = (v1 * v2) / (v1.Length() * v2.Length());
  if (co > 1) co = 1;
  if (co < -1) co = -1;
  return acos ( co );
}

void Vec3d :: GetNormal (Vec3d & n) const
  {
  if (fabs (X()) > fabs (Z()))
    {
    n.X() = -Y();
    n.Y() = X();
    n.Z() = 0;
    }
  else
    {
    n.X() = 0;
    n.Y() = Z();
    n.Z() = -Y();
    }
  double len = n.Length();
  if (len == 0)
    {
    n.X() = 1;
    n.Y() = n.Z() = 0;
    }
  else
    n /= len;
  }

Box3d :: Box3d ( double aminx, double amaxx,
                 double aminy, double amaxy,
                 double aminz, double amaxz )
{
  minx[0] = aminx; maxx[0] = amaxx;
  minx[1] = aminy; maxx[1] = amaxy;
  minx[2] = aminz; maxx[2] = amaxz;
}

Box3d :: Box3d ( const Box3d & b2 )
{
  for (int i = 0; i < 3; i++)
    {
      minx[i] = b2.minx[i];
      maxx[i] = b2.maxx[i];
    }
}

Box3d :: Box3d ( const Box<3> & b2 )
{
  for (int i = 0; i < 3; i++)
    {
      minx[i] = b2.PMin()(i);
      maxx[i] = b2.PMax()(i);
    }
}

void Box3d :: GetPointNr (int i, Point3d & point) const
{
  i--;
  point.X() = (i & 1) ? maxx[0] : minx[0];
  point.Y() = (i & 2) ? maxx[1] : minx[1];
  point.Z() = (i & 4) ? maxx[2] : minx[2];
}


void Box3d :: Increase (double d)
{
  for (int i = 0; i <= 2; i++)
    {
      minx[i] -= d;
      maxx[i] += d;
    }
}

void Box3d :: IncreaseRel (double /* rel */)
{
  for (int i = 0; i <= 2; i++)
    {
      double d = 0.5 * (maxx[i] - minx[i]);
      minx[i] -= d;
      maxx[i] += d;
    }
}

Box3d :: Box3d (const Point3d& p1, const Point3d& p2)
{
  minx[0] = min2 (p1.X(), p2.X());
  minx[1] = min2 (p1.Y(), p2.Y());
  minx[2] = min2 (p1.Z(), p2.Z());
  maxx[0] = max2 (p1.X(), p2.X());
  maxx[1] = max2 (p1.Y(), p2.Y());
  maxx[2] = max2 (p1.Z(), p2.Z());
}

const Box3d& Box3d :: operator+=(const Box3d& b)
{
  minx[0] = min2 (minx[0], b.minx[0]);
  minx[1] = min2 (minx[1], b.minx[1]);
  minx[2] = min2 (minx[2], b.minx[2]);
  maxx[0] = max2 (maxx[0], b.maxx[0]);
  maxx[1] = max2 (maxx[1], b.maxx[1]);
  maxx[2] = max2 (maxx[2], b.maxx[2]);

  return *this;
}

Point3d Box3d :: MaxCoords() const
{
  return Point3d(maxx[0], maxx[1], maxx[2]);
}

Point3d Box3d :: MinCoords() const
{
  return Point3d(minx[0], minx[1], minx[2]);
}

void Box3d :: WriteData(ofstream& fout) const
{
  for(int i = 0; i < 3; i++)
    {
      fout << minx[i] << " " << maxx[i] << " ";
    }
  fout << "\n";
}

void Box3d :: ReadData(ifstream& fin)
{
  for(int i = 0; i < 3; i++)
    {
      fin >> minx[i];
      fin >> maxx[i];
    }
}

Box3dSphere :: Box3dSphere ( double aminx, double amaxx,
			     double aminy, double amaxy,
			     double aminz, double amaxz )
  : Box3d (aminx, amaxx, aminy, amaxy, aminz, amaxz)
{
  CalcDiamCenter ();
}

void Box3dSphere :: CalcDiamCenter ()
{
  diam = sqrt( sqr (maxx[0] - minx[0]) +
	       sqr (maxx[1] - minx[1]) + 
	       sqr (maxx[2] - minx[2]));
  
  c.X() = 0.5 * (minx[0] + maxx[0]);
  c.Y() = 0.5 * (minx[1] + maxx[1]);
  c.Z() = 0.5 * (minx[2] + maxx[2]);
  
  inner = min2 ( min2 (maxx[0] - minx[0], maxx[1] - minx[1]), maxx[2] - minx[2]) / 2;
}

void Box3dSphere :: GetSubBox (int i, Box3dSphere & sbox) const
{
  i--;
  if (i & 1)
    {
      sbox.minx[0] = c.X();
      sbox.maxx[0] = maxx[0];
    }
  else
    {
      sbox.minx[0] = minx[0];
      sbox.maxx[0] = c.X();
    }
  if (i & 2)
    {
      sbox.minx[1] = c.Y();
      sbox.maxx[1] = maxx[1];
    }
  else
    {
      sbox.minx[1] = minx[1];
      sbox.maxx[1] = c.Y();
    }
  if (i & 4)
    {
      sbox.minx[2] = c.Z();
      sbox.maxx[2] = maxx[2];
    }
  else
    {
      sbox.minx[2] = minx[2];
      sbox.maxx[2] = c.Z();
    }
  
  sbox.c.X() = 0.5 * (sbox.minx[0] + sbox.maxx[0]);
  sbox.c.Y() = 0.5 * (sbox.minx[1] + sbox.maxx[1]);
  sbox.c.Z() = 0.5 * (sbox.minx[2] + sbox.maxx[2]);
  sbox.diam = 0.5 * diam;
  sbox.inner = 0.5 * inner;
}

void Transpose (Vec3d & v1, Vec3d & v2, Vec3d & v3)
{
  Swap (v1.Y(), v2.X());
  Swap (v1.Z(), v3.X());
  Swap (v2.Z(), v3.Y());
}

int SolveLinearSystem (const Vec3d & col1, const Vec3d & col2,
		       const Vec3d & col3, const Vec3d & rhs,
		       Vec3d & sol)
{
  // changed by MW
  double matrix[3][3];
  double locrhs[3];
  int retval = 0;

  for(int i=0; i<3; i++)
    {
      matrix[i][0] = col1.X(i+1);
      matrix[i][1] = col2.X(i+1);
      matrix[i][2] = col3.X(i+1);
      locrhs[i] = rhs.X(i+1);
    }

  for(int i=0; i<2; i++)
    {
      int pivot = i;
      double maxv = fabs(matrix[i][i]);
      for(int j=i+1; j<3; j++)
	if(fabs(matrix[j][i]) > maxv)
	  {
	    maxv = fabs(matrix[j][i]);
	    pivot = j;
	  }

      if(fabs(maxv) > 1e-40)
	{
	  if(pivot != i)
	    {
	      swap(matrix[i][0],matrix[pivot][0]);
	      swap(matrix[i][1],matrix[pivot][1]);
	      swap(matrix[i][2],matrix[pivot][2]);
	      swap(locrhs[i],locrhs[pivot]);
	    }
	  for(int j=i+1; j<3; j++)
	    {
	      double fac = matrix[j][i] / matrix[i][i];
	      
	      for(int k=i+1; k<3; k++)
		matrix[j][k] -= fac*matrix[i][k];
	      locrhs[j] -= fac*locrhs[i];
	    }
	}
      else
	retval = 1;
    }

  if(fabs(matrix[2][2]) < 1e-40)
    retval = 1;

  if(retval != 0)
    return retval;
  

  for(int i=2; i>=0; i--)
    {
      double sum = locrhs[i];
      for(int j=2; j>i; j--)
	sum -= matrix[i][j]*sol.X(j+1);

      sol.X(i+1) = sum/matrix[i][i];
    }

  return 0;
  
}

int PseudoInverse (const Vec3d & col1,
		   const Vec3d & col2,
		   Vec3d & inv1,
		   Vec3d & inv2)
{
  double a11 = col1 * col1;
  double a12 = col1 * col2;
  double a22 = col2 * col2;
  
  double det = a11 * a22 - a12 * a12;

  if (fabs (det) < 1e-12 * col1.Length() * col2.Length())
    {
      inv1 = Vec3d (0, 0, 0);
      inv2 = Vec3d (0, 0, 0);
      return 1;
    }

  double ia11 = a22 / det;
  double ia12 = -a12 / det;
  double ia22 = a11 / det;

  inv1 = ia11 * col1 + ia12 * col2;
  inv2 = ia12 * col1 + ia22 * col2;

  return 0;
}

QuadraticFunction3d :: 
QuadraticFunction3d (const Point3d & p, const Vec3d & v)
{
  Vec3d hv(v);
  hv /= (hv.Length() + 1e-12);
  Vec3d t1, t2;
  hv.GetNormal (t1);
  Cross (hv, t1, t2);
  
  double t1p = t1.X() * p.X() + t1.Y() * p.Y() + t1.Z() * p.Z();
  double t2p = t2.X() * p.X() + t2.Y() * p.Y() + t2.Z() * p.Z();
  c0 = sqr (t1p) + sqr (t2p);
  cx = -2 * (t1p * t1.X() + t2p * t2.X());
  cy = -2 * (t1p * t1.Y() + t2p * t2.Y());
  cz = -2 * (t1p * t1.Z() + t2p * t2.Z());

  cxx = t1.X() * t1.X() + t2.X() * t2.X();
  cyy = t1.Y() * t1.Y() + t2.Y() * t2.Y();
  czz = t1.Z() * t1.Z() + t2.Z() * t2.Z();

  cxy = 2 * t1.X() * t1.Y() + 2 * t2.X() * t2.Y();
  cxz = 2 * t1.X() * t1.Z() + 2 * t2.X() * t2.Z();
  cyz = 2 * t1.Y() * t1.Z() + 2 * t2.Y() * t2.Z();
}

void referencetransform :: Set (const Point3d & p1, const Point3d & p2,
                                const Point3d & p3, double ah)
{
  ex = p2 - p1;
  ex /= ex.Length();
  ey = p3 - p1;
  ey -= (ex * ey) * ex;
  ey /= ey.Length();
  ez = Cross (ex, ey);
  rp = p1;
  h = ah;

  exh = ah * ex;
  eyh = ah * ey;
  ezh = ah * ez;
  ah = 1 / ah;
  ex_h = ah * ex;
  ey_h = ah * ey;
  ez_h = ah * ez;
}

void referencetransform :: ToPlain (const Point3d & p, Point3d & pp) const
{
  Vec3d v;
  v = p - rp;
  pp.X() = (ex_h * v);
  pp.Y() = (ey_h * v);
  pp.Z() = (ez_h * v);
}

void referencetransform :: ToPlain (const Array<Point3d> & p,
                                    Array<Point3d> & pp) const
{
  Vec3d v;
  int i;

  pp.SetSize (p.Size());
  for (i = 1; i <= p.Size(); i++)
    {
      v = p.Get(i) - rp;
      pp.Elem(i).X() = (ex_h * v);
      pp.Elem(i).Y() = (ey_h * v);
      pp.Elem(i).Z() = (ez_h * v);
    }
}

void referencetransform :: FromPlain (const Point3d & pp, Point3d & p) const
{
  Vec3d v;
  v.X() = pp.X() * exh.X() + pp.Y() * eyh.X() + pp.Z() * ezh.X();
  v.Y() = pp.X() * exh.Y() + pp.Y() * eyh.Y() + pp.Z() * ezh.Y();
  v.Z() = pp.X() * exh.Z() + pp.Y() * eyh.Z() + pp.Z() * ezh.Z();
  p.X() = rp.X() + v.X();
  p.Y() = rp.Y() + v.Y();
  p.Z() = rp.Z() + v.Z();
}
}