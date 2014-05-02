#pragma once
/* *************************************************************************/
/* File:   geomtest3d.hh                                                   */
/* Author: Joachim Schoeberl                                               */
/* Date:   13. Feb. 98                                                     */
/* *************************************************************************/


namespace netgen
{

extern void
LocalCoordinates (const Vec3d & e1, const Vec3d & e2,
		  const Vec3d & v, double & lam1, double & lam2);

/// Minimal distance of point p to the line segment [lp1,lp2]
extern double MinDistLP2 (const Point2d & lp1, const Point2d & lp2, const Point2d & p);

/// Minimal distance of point p to the line segment [lp1,lp2]
extern double MinDistLP2 (const Point3d & lp1, const Point3d & lp2, const Point3d & p);

/// Minimal distance of point p to the triangle segment [tp1,tp2,pt3]
extern double MinDistTP2 (const Point3d & tp1, const Point3d & tp2, 
			  const Point3d & tp3, const Point3d & p);

}