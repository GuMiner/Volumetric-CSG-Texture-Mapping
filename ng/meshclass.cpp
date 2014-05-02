#include "stdafx.h" 
#include "mystdlib.h"
#include "meshing.hpp"

namespace netgen
{
    // Meshing parameters for global access
    extern MeshingParameters mparam;

  Mesh :: Mesh ()
    : surfarea(*this)
  {
    boundaryedges = NULL;
    surfelementht = NULL; 
    segmentht = NULL;

    lochfunc = NULL;
    hglob = 1e10;
    hmin = 0;
    numvertices = -1;

    ident = new Identifications (*this);

    coarsemesh = NULL;

    geomtype = NO_GEOM;

    bcnames.SetSize(0);
  }

  Mesh :: ~Mesh()
  {
    delete lochfunc;
    delete boundaryedges;
    delete surfelementht;
    delete segmentht;
    delete ident;
    delete coarsemesh;

    for (int i = 0; i < bcnames.Size(); i++ )
      if ( bcnames[i] ) delete bcnames[i];

  }

  Mesh & Mesh :: operator= (const Mesh & mesh2)
  {
    points = mesh2.points;
    // eltyps = mesh2.eltyps;
    segments = mesh2.segments;
    surfelements = mesh2.surfelements;
    lockedpoints = mesh2.lockedpoints;
    facedecoding = mesh2.facedecoding;

    bcnames.SetSize( mesh2.bcnames.Size() );
    for ( int i = 0; i < mesh2.bcnames.Size(); i++ )
      if ( mesh2.bcnames[i] ) bcnames[i] = new string ( *mesh2.bcnames[i] );
      else bcnames[i] = 0;

    return *this;
  }

  void Mesh :: DeleteMesh()
  {
    points.SetSize(0);
    segments.SetSize(0);
    surfelements.SetSize(0);
    lockedpoints.SetSize(0);
    surfacesonnode.SetSize(0);

    delete boundaryedges;
    boundaryedges = NULL;

    openelements.SetSize(0);
    facedecoding.SetSize(0);

    delete ident;
    ident = new Identifications (*this);

    for ( int i = 0; i < bcnames.Size(); i++ )
      if ( bcnames[i] ) delete bcnames[i];
  }

  void Mesh :: ClearSurfaceElements()
  { 
    surfelements.SetSize(0); 
    for (int i = 0; i < facedecoding.Size(); i++)
      facedecoding[i].firstelement = -1;
  }

  PointIndex Mesh :: AddPoint (const Point3d & p, int layer)
  { 
    return AddPoint (p, layer, INNERPOINT);
  }

  PointIndex Mesh :: AddPoint (const Point3d & p, int layer, POINTTYPE type)
  { 
    PointIndex pi = points.End();
    points.Append ( MeshPoint (p, layer, type) ); 

    return pi;
  }

  SegmentIndex Mesh :: AddSegment (const Segment & s)
  { 
    int maxn = max2 (s[0], s[1]);
    maxn += 1-PointIndex::BASE;

    if (maxn <= points.Size())
      {
        if (points[s[0]].Type() > EDGEPOINT)
          points[s[0]].SetType (EDGEPOINT);
        if (points[s[1]].Type() > EDGEPOINT)
          points[s[1]].SetType (EDGEPOINT);
      }

    SegmentIndex si = segments.Size();
    segments.Append (s); 

    return si;
  }

  SurfaceElementIndex Mesh :: AddSurfaceElement (const Element2d & el)
  {
    int maxn = el[0];
    for (int i = 1; i < el.GetNP(); i++)
      if (el[i] > maxn) maxn = el[i];

    maxn += 1-PointIndex::BASE;

    if (maxn <= points.Size())
      {
        for (int i = 0; i < el.GetNP(); i++)
          if (points[el[i]].Type() > SURFACEPOINT)
            points[el[i]].SetType(SURFACEPOINT);
      }

    SurfaceElementIndex si = surfelements.Size();
    surfelements.Append (el); 

    if (el.index > facedecoding.Size())
      cerr << "has no facedecoding: fd.size = " << facedecoding.Size() << ", ind = " << el.index << endl;

    surfelements.Last().next = facedecoding[el.index-1].firstelement;
    facedecoding[el.index-1].firstelement = si;

    if (SurfaceArea().Valid())
      SurfaceArea().Add (el);

    return si;
  }


  void Mesh :: SetAllocSize(int nnodes, int nsegs, int nsel, int nel)
  {
    points.SetAllocSize(nnodes);
    segments.SetAllocSize(nsegs);
    surfelements.SetAllocSize(nsel);
  }

  void Mesh :: BuildBoundaryEdges(void)
  {
    delete boundaryedges;

    boundaryedges = new INDEX_2_CLOSED_HASHTABLE<int>
      (3 * (GetNSE() + GetNOpenElements()) + GetNSeg() + 1);

    for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
      {
        const Element2d & sel = surfelements[sei];
        if (sel.IsDeleted()) continue;

        // int si = sel.GetIndex();

        for (int j = 0; j < sel.GetNP(); j++)
          {
            INDEX_2 i2;
            i2.I1() = sel.PNumMod(j+1);
            i2.I2() = sel.PNumMod(j+2);
            i2.Sort();
            if (sel.GetNP() <= 4)
              boundaryedges->Set (i2, 1);
          }
      }

    for (int i = 0; i < openelements.Size(); i++)
      {
        const Element2d & sel = openelements[i];
        for (int j = 0; j < sel.GetNP(); j++)
          {
            INDEX_2 i2;
            i2.I1() = sel.PNumMod(j+1);
            i2.I2() = sel.PNumMod(j+2);
            i2.Sort();
            boundaryedges->Set (i2, 1);

            points[sel[j]].SetType(FIXEDPOINT);
          }
      }

    for (int i = 0; i < GetNSeg(); i++)
      {
        const Segment & seg = segments[i];
        INDEX_2 i2(seg[0], seg[1]);
        i2.Sort();

        boundaryedges -> Set (i2, 2);
      }
  }

  void Mesh :: CalcSurfacesOfNode ()
  {
    surfacesonnode.SetSize (GetNP());

    delete boundaryedges;
    boundaryedges = NULL;

    delete surfelementht;
    delete segmentht;

    surfelementht = new INDEX_3_CLOSED_HASHTABLE<int> (3*GetNSE() + 1);
    segmentht = new INDEX_2_CLOSED_HASHTABLE<int> (3*GetNSeg() + 1);

    for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
      {
        const Element2d & sel = surfelements[sei];
        if (sel.IsDeleted()) continue;

        int si = sel.GetIndex();

        for (int j = 0; j < sel.GetNP(); j++)
          {
            PointIndex pi = sel[j];
            bool found = 0;
            for (int k = 0; k < surfacesonnode[pi].Size(); k++)
              if (surfacesonnode[pi][k] == si)
                {
                  found = 1;
                  break;
                }

            if (!found)
              surfacesonnode.Add (pi, si);
          }
      }
    for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
      {
        const Element2d & sel = surfelements[sei];
        if (sel.IsDeleted()) continue;

        INDEX_3 i3;
        i3.I1() = sel.PNum(1);
        i3.I2() = sel.PNum(2);
        i3.I3() = sel.PNum(3);
        i3.Sort();
        surfelementht -> Set (i3, sei);   // war das wichtig ???    sel.GetIndex());
      }

    int np = GetNP();

    for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
        points[pi].SetType (INNERPOINT);

    if (GetNFD() == 0) 
        {
        for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
            {
            const Element2d & sel = surfelements[sei];
            if (sel.IsDeleted()) continue;
            for (int j = 0;  j < sel.GetNP(); j++)
                {
                PointIndex pi = SurfaceElement(sei)[j];
                points[pi].SetType(FIXEDPOINT);
                }
            }
        }
    else
        {
        for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
            {
            const Element2d & sel = surfelements[sei];
            if (sel.IsDeleted()) continue;
            for (int j = 0; j < sel.GetNP(); j++)
                {
                PointIndex pi = sel[j];
                int ns = surfacesonnode[pi].Size();
                if (ns == 1)
                    points[pi].SetType(SURFACEPOINT);
                if (ns == 2)
                    points[pi].SetType(EDGEPOINT);
                if (ns >= 3)
                    points[pi].SetType(FIXEDPOINT);
                }      
            }
        }

    for (int i = 0; i < segments.Size(); i++)
      {
        const Segment & seg = segments[i];
        for (int j = 1; j <= 2; j++)
          {
            PointIndex hi = (j == 1) ? seg[0] : seg[1];
            if (points[hi].Type() == INNERPOINT ||
                points[hi].Type() == SURFACEPOINT)
              points[hi].SetType(EDGEPOINT);
          }
      }
    
    for (int i = 0; i < lockedpoints.Size(); i++)
      points[lockedpoints[i]].SetType(FIXEDPOINT);

    for (int i = 0; i < GetNSeg(); i++)
      {
        const Segment & seg = segments[i];
        INDEX_2 i2(seg[0], seg[1]);
        i2.Sort();

        segmentht -> Set (i2, i);
      }
  }

  void Mesh :: FixPoints (const BitArray & fixpoints)
  {
    if (fixpoints.Size() != GetNP())
      {
        cerr << "Mesh::FixPoints: sizes don't fit" << endl;
        return;
      }
    int np = GetNP();
    for (int i = 1; i <= np; i++)
      if (fixpoints.Test(i))
        {
          points.Elem(i).SetType (FIXEDPOINT);
        }
  }

  void Mesh :: SetLocalH (const Point3d & pmin, const Point3d & pmax, double grading)
  {
    Point3d c = Center (pmin, pmax);
    double d = max3 (pmax.X()-pmin.X(),
                     pmax.Y()-pmin.Y(),
                     pmax.Z()-pmin.Z());
    d /= 2;
    Point3d pmin2 = c - Vec3d (d, d, d);
    Point3d pmax2 = c + Vec3d (d, d, d);


    delete lochfunc;
    lochfunc = new LocalH (pmin2, pmax2, grading);
  }

  void Mesh :: RestrictLocalH (const Point3d & p, double hloc)
  {
    if(hloc < hmin)
      hloc = hmin;

    if (!lochfunc)
      {
        Point3d boxmin, boxmax;
        GetBox (boxmin, boxmax);
        SetLocalH (boxmin, boxmax, 0.8f);
      }

    lochfunc -> SetH (p, hloc);
  }

  void Mesh :: RestrictLocalHLine (const Point3d & p1, 
                                   const Point3d & p2,
                                   double hloc)
  {
    if(hloc < hmin)
      hloc = hmin;

    int i;
    int steps = int (Dist (p1, p2) / hloc) + 2;
    Vec3d v(p1, p2);

    for (i = 0; i <= steps; i++)
      {
        Point3d p = p1 + (double(i)/double(steps) * v);
        RestrictLocalH (p, hloc);
      }
  }

  void Mesh :: SetMinimalH (double h)
  {
    hmin = h;
  }

  void Mesh :: SetGlobalH (double h)
  {
    hglob = h;
  }

  double Mesh :: MaxHDomain (int dom) const
  {
    if (maxhdomain.Size())
      return maxhdomain.Get(dom);
    else
      return 1e10;
  }

  void Mesh :: SetMaxHDomain (const Array<double> & mhd)
  {
    maxhdomain.SetSize(mhd.Size());
    for (int i = 1; i <= mhd.Size(); i++)
      maxhdomain.Elem(i) = mhd.Get(i);
  }

  double Mesh :: GetH (const Point3d & p) const
  {
    double hmin = hglob;
    if (lochfunc)
      {
        double hl = lochfunc->GetH (p);
        if (hl < hglob)
          hmin = hl;
      }
    return hmin;
  }

  double Mesh :: GetMinH (const Point3d & pmin, const Point3d & pmax)
  {
    double hmin = hglob;
    if (lochfunc)
      {
        double hl = lochfunc->GetMinH (pmin, pmax);
        if (hl < hmin)
          hmin = hl;
      }
    return hmin;
  }

  double Mesh :: AverageH (int surfnr) const
  {
    int i, j, n;
    double hi, hsum;
    double maxh = 0, minh = 1e10;

    hsum = 0;
    n = 0;
    for (i = 1; i <= GetNSE(); i++)
      {
        const Element2d & el = SurfaceElement(i);
        if (surfnr == 0 || el.GetIndex() == surfnr)
          {
            for (j = 1; j <= 3; j++)
              {
                hi = Dist (Point (el.PNumMod(j)), 
                           Point (el.PNumMod(j+1)));

                hsum += hi;

                if (hi > maxh) maxh = hi;
                if (hi < minh) minh = hi;
                n++;
              }
          }
      }

    cout <<  "minh = " << minh <<  " avh = "<< (hsum/n)<< " maxh = "<< maxh<<endl;
    return (hsum / n);
  }

  void Mesh :: CalcLocalH (double grading) 
  {
    if (!lochfunc)
      {
        Point3d pmin, pmax;
        GetBox (pmin, pmax);
	SetLocalH (pmin, pmax, grading);
      }

    for (int i = 0; i < GetNSE(); i++)
      {
        const Element2d & el = surfelements[i];
        int j;

        if (el.GetNP() == 3)
          {
            double hel = -1;
            for (j = 1; j <= 3; j++)
              {
                const Point3d & p1 = points[el.PNumMod(j)];
                const Point3d & p2 = points[el.PNumMod(j+1)];

                if (!ident -> UsedSymmetric (el.PNumMod(j),
                                             el.PNumMod(j+1)))
                  {
                    double hedge = Dist (p1, p2);
                    if (hedge > hel)
                      hel = hedge;
                  }
              }

            if (hel > 0)
              {
                const Point3d & p1 = points[el.PNum(1)];
                const Point3d & p2 = points[el.PNum(2)];
                const Point3d & p3 = points[el.PNum(3)];
                lochfunc->SetH (Center (p1, p2, p3), hel);
              }
          }
        else
          {
            {
              const Point3d & p1 = points[el.PNum(1)];
              const Point3d & p2 = points[el.PNum(2)];
              lochfunc->SetH (Center (p1, p2), 2 * Dist (p1, p2));
            }
            {
              const Point3d & p1 = points[el.PNum(3)];
              const Point3d & p2 = points[el.PNum(4)];
              lochfunc->SetH (Center (p1, p2), 2 * Dist (p1, p2));
            }
          }
      }

    for (int i = 0; i < GetNSeg(); i++)
      {
        const Segment & seg = segments[i];
        const Point3d & p1 = points[seg[0]];
        const Point3d & p2 = points[seg[1]];
        if (!ident -> UsedSymmetric (seg[0], seg[1]))
          {
            lochfunc->SetH (Center (p1, p2), Dist (p1, p2));
          }
      }
  }

  void Mesh :: RestrictLocalH (resthtype rht, int nr, double loch)
  {
    int i;
    switch (rht)
      {
      case RESTRICTH_FACE:
        {
          for (i = 1; i <= GetNSE(); i++)
            {
              const Element2d & sel = SurfaceElement(i);
              if (sel.GetIndex() == nr)
                RestrictLocalH (RESTRICTH_SURFACEELEMENT, i, loch);
            }
          break;
        }
      case RESTRICTH_EDGE:
        {
          for (i = 1; i <= GetNSeg(); i++)
            {
              const Segment & seg = LineSegment(i);
              if (seg.edgenr == nr)
                RestrictLocalH (RESTRICTH_SEGMENT, i, loch);
            }
          break;
        }
      case RESTRICTH_POINT:
        {
          RestrictLocalH (Point (nr), loch);
          break;
        }

      case RESTRICTH_SURFACEELEMENT:
        {
          const Element2d & sel = SurfaceElement(nr);
          Point3d p = Center (Point(sel.PNum(1)),
                              Point(sel.PNum(2)),
                              Point(sel.PNum(3)));
          RestrictLocalH (p, loch);
          break;
        }
      case RESTRICTH_SEGMENT:
        {
          const Segment & seg = LineSegment(nr);
          RestrictLocalHLine (Point (seg[0]), Point(seg[1]), loch);
          break;
        }
      }
  }

  void Mesh :: GetBox (Point3d & pmin, Point3d & pmax, int dom) const
  {
    if (points.Size() == 0)
      {
        pmin = pmax = Point3d(0,0,0);
        return;
      }

    if (dom <= 0)
      {
        pmin = Point3d (1e10, 1e10, 1e10);
        pmax = Point3d (-1e10, -1e10, -1e10); 

        for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
          {
            pmin.SetToMin ( (*this) [pi] );
            pmax.SetToMax ( (*this) [pi] );
          }
      }
    else
      {
        int j, nse = GetNSE();
        SurfaceElementIndex sei;

        pmin = Point3d (1e10, 1e10, 1e10);
        pmax = Point3d (-1e10, -1e10, -1e10); 
        for (sei = 0; sei < nse; sei++)
          {
            const Element2d & el = (*this)[sei];
            if (el.IsDeleted() ) continue;

            if (dom == -1 || el.GetIndex() == dom)
              {
                for (j = 0; j < 3; j++)
                  {
                    pmin.SetToMin ( (*this) [el[j]] );
                    pmax.SetToMax ( (*this) [el[j]] );
                  }
              }
          }
      }

    if (pmin.X() > 0.5e10)
      {
        pmin = pmax = Point3d(0,0,0);
      }
  }

  void Mesh :: GetBox (Point3d & pmin, Point3d & pmax, POINTTYPE ptyp) const
  {
    if (points.Size() == 0)
      {
        pmin = pmax = Point3d(0,0,0);
        return;
      }

    pmin = Point3d (1e10, 1e10, 1e10);
    pmax = Point3d (-1e10, -1e10, -1e10); 

    for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
      if (points[pi].Type() <= ptyp)
        {
          pmin.SetToMin ( (*this) [pi] );
          pmax.SetToMax ( (*this) [pi] );
        }
  }

  void Mesh :: AddLockedPoint (PointIndex pi)
  { 
    lockedpoints.Append (pi); 
  }

  void Mesh :: ClearLockedPoints ()
  { 
    lockedpoints.SetSize (0); 
  }

  void Mesh :: Compress ()
  {
    Array<PointIndex,PointIndex::BASE,PointIndex> op2np(GetNP());
    Array<MeshPoint> hpoints;
    BitArrayChar<PointIndex::BASE> pused(GetNP());

    for (int i = 0; i < surfelements.Size(); i++)
      if (surfelements[i].IsDeleted())
        {
          surfelements.Delete(i);
          i--;
        }

    for (int i = 0; i < segments.Size(); i++)
      if (segments[i][0] <= PointIndex::BASE-1)
        {
          segments.Delete(i);
          i--;
        }

    pused.Clear();

    for (int i = 0; i < surfelements.Size(); i++)
      {
        const Element2d & el = surfelements[i];
        for (int j = 0; j < el.GetNP(); j++)
          pused.Set (el[j]);
      }

    for (int i = 0; i < segments.Size(); i++)
      {
        const Segment & seg = segments[i];
        pused.Set (seg[0]);
        pused.Set (seg[1]);
      }

    for (int i = 0; i < openelements.Size(); i++)
      {
        const Element2d & el = openelements[i];
        for (int j = 0; j < el.GetNP(); j++)
          pused.Set(el[j]);
      }

    for (int i = 0; i < lockedpoints.Size(); i++)
      pused.Set (lockedpoints[i]);

    int npi = PointIndex::BASE-1;

    for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
      if (pused.Test(pi))
        {
          npi++;
          op2np[pi] = npi;
          hpoints.Append (points[pi]);
        }
      else
        op2np[pi] = -1;

    points.SetSize(0);
    for (int i = 0; i < hpoints.Size(); i++)
      points.Append (hpoints[i]);

    for (int i = 1; i <= surfelements.Size(); i++)
      {
        Element2d & el = SurfaceElement(i);
        for (int j = 0; j < el.GetNP(); j++)
          el[j] = op2np[el[j]];
      }

    for (int i = 0; i < segments.Size(); i++)
      {
        Segment & seg = segments[i];
        seg[0] = op2np[seg[0]];
        seg[1] = op2np[seg[1]];
      }

    for (int i = 1; i <= openelements.Size(); i++)
      {
        Element2d & el = openelements.Elem(i);
        for (int j = 0; j < el.GetNP(); j++)
          el[j] = op2np[el[j]];
      }  


    for (int i = 0; i < lockedpoints.Size(); i++)
      lockedpoints[i] = op2np[lockedpoints[i]];

    for (int i = 0; i < facedecoding.Size(); i++)
      facedecoding[i].firstelement = -1;
    for (int i = surfelements.Size()-1; i >= 0; i--)
      {
        int ind = surfelements[i].GetIndex();
        surfelements[i].next = facedecoding[ind-1].firstelement;
        facedecoding[ind-1].firstelement = i;
      }

    CalcSurfacesOfNode();
  }
  
  bool Mesh :: LegalTrig (const Element2d & el) const
  {
    return 1;
    if ( /* hp */ 1)  // needed for old, simple hp-refinement
      { 
        // trigs with 2 or more segments are illegal
        int i;
        int nseg = 0;

        if (!segmentht)
          {
            cerr << "no segmentht allocated" << endl;
            return 0;
          }

        //      Point3d cp(0.5, 0.5, 0.5);
        for (i = 1; i <= 3; i++)
          {
            INDEX_2 i2(el.PNumMod (i), el.PNumMod (i+1));
            i2.Sort();
            if (segmentht -> Used (i2))
              nseg++;
          }
        if (nseg >= 2) 
          return 0;
      }
    return 1;
  }

  void Mesh :: SurfaceMeshOrientation ()
  {
    int i, j;
    int nse = GetNSE();

    BitArray used(nse);
    used.Clear();
    INDEX_2_HASHTABLE<int> edges(nse+1);

    bool haschanged = 0;

    const Element2d & tri = SurfaceElement(1);
    for (j = 1; j <= 3; j++)
      {
        INDEX_2 i2(tri.PNumMod(j), tri.PNumMod(j+1));
        edges.Set (i2, 1);
      }
    used.Set(1);

    bool unused;
    do
      {
        bool changed;
        do
          {
            changed = 0;
            for (i = 1; i <= nse; i++)
              if (!used.Test(i))
                {
                  Element2d & el = surfelements.Elem(i);
                  int found = 0, foundrev = 0;
                  for (j = 1; j <= 3; j++)
                    {
                      INDEX_2 i2(el.PNumMod(j), el.PNumMod(j+1));
                      if (edges.Used(i2))
                        foundrev = 1;
                      swap (i2.I1(), i2.I2());
                      if (edges.Used(i2))
                        found = 1;
                    }

                  if (found || foundrev)
                    {
                      if (foundrev)
                        swap (el.PNum(2), el.PNum(3));

                      changed = 1;
                      for (j = 1; j <= 3; j++)
                        {
                          INDEX_2 i2(el.PNumMod(j), el.PNumMod(j+1));
                          edges.Set (i2, 1);
                        }
                      used.Set (i);
                    }
                }
            if (changed)
              haschanged = 1;
          }
        while (changed);

        unused = 0;
        for (i = 1; i <= nse; i++)
          if (!used.Test(i))
            {
              unused = 1;
              const Element2d & tri = SurfaceElement(i);
              for (j = 1; j <= 3; j++)
                {
                  INDEX_2 i2(tri.PNumMod(j), tri.PNumMod(j+1));
                  edges.Set (i2, 1);
                }
              used.Set(i);
              break;
            }
      }
    while (unused);
  }

  void Mesh :: SplitSeparatedFaces ()
  {
      if (mparam.enableOutput)
        cout << "Split Separate Faces" << endl;
    int fdi;
    int np = GetNP();

    BitArray usedp(np);
    Array<SurfaceElementIndex> els_of_face;

    fdi = 1;
    while (fdi <= GetNFD())
      {
        GetSurfaceElementsOfFace (fdi, els_of_face);

        if (els_of_face.Size() == 0) continue;

        SurfaceElementIndex firstel = els_of_face[0];

        usedp.Clear();
        for (int j = 1; j <= SurfaceElement(firstel).GetNP(); j++)
          usedp.Set (SurfaceElement(firstel).PNum(j));

        bool changed;
        do
          {
            changed = false;

            for (int i = 0; i < els_of_face.Size(); i++)
              {
                const Element2d & el = SurfaceElement(els_of_face[i]);

                bool has = 0;
                bool hasno = 0;
                for (int j = 0; j < el.GetNP(); j++)
                  {
                    if (usedp.Test(el[j]))
                      has = true;
                    else
                      hasno = true;
                  }

                if (has && hasno)
                  changed = true;

                if (has)
                  for (int j = 0; j < el.GetNP(); j++)
                    usedp.Set (el[j]);
              }
          }
        while (changed);

        int nface = 0;
        for (int i = 0; i < els_of_face.Size(); i++)
          {
            Element2d & el = SurfaceElement(els_of_face[i]);

            int hasno = 0;
            for (int j = 1; j <= el.GetNP(); j++)
              if (!usedp.Test(el.PNum(j)))
                hasno = 1;

            if (hasno)
              {
                if (!nface)
                  {
                    FaceDescriptor nfd = GetFaceDescriptor(fdi);
                    nface = AddFaceDescriptor (nfd);
                  }

                el.SetIndex (nface);
              }
          }

        // reconnect list
        if (nface)
          {
            facedecoding[nface-1].firstelement = -1;
            facedecoding[fdi-1].firstelement = -1;

            for (int i = 0; i < els_of_face.Size(); i++)
              {
                int ind = SurfaceElement(els_of_face[i]).GetIndex();
                SurfaceElement(els_of_face[i]).next = facedecoding[ind-1].firstelement;
                facedecoding[ind-1].firstelement = els_of_face[i];
              }
          }

        fdi++;
      }
  }

  void Mesh :: RebuildSurfaceElementLists ()
  {
    for (int i = 0; i < facedecoding.Size(); i++)
      facedecoding[i].firstelement = -1;
    for (int i = surfelements.Size()-1; i >= 0; i--)
      {
        int ind = surfelements[i].GetIndex();
        surfelements[i].next = facedecoding[ind-1].firstelement;
        facedecoding[ind-1].firstelement = i;
      }
  }

  void Mesh :: GetSurfaceElementsOfFace (int facenr, Array<SurfaceElementIndex> & sei) const
  {
     /* Philippose - 01/10/2009
     Commented out the following lines, and activated the originally 
     commented out lines above because of a bug which causes corruption 
     of the variable "facedecoding" when a mesh is converted to second order
     */

     sei.SetSize(0);

     SurfaceElementIndex si = facedecoding[facenr-1].firstelement;
     while (si != -1)
     {
        if ( (*this)[si].GetIndex () == facenr && (*this)[si][0] >= PointIndex::BASE &&
             !(*this)[si].IsDeleted() )
        {
           sei.Append (si);
        }

        si = (*this)[si].next;
     }
  }

  void Mesh :: ComputeNVertices ()
  {
    numvertices = 0;
    for (int i = 1; i <= GetNSE(); i++)
      {
        const Element2d & el = SurfaceElement(i);
        for (int j = 1; j <= el.GetNV(); j++)
          if (el.PNum(j) > numvertices)
            numvertices = el.PNum(j);
      } 

    numvertices += 1- PointIndex::BASE;
  }

  int Mesh :: GetNV () const
  {
    if (numvertices < 0)
      return GetNP();
    else
      return numvertices;
  }

  bool Mesh :: PureTrigMesh (int faceindex) const
  {
    if (!faceindex)
      {
	for (int i = 1; i <= GetNSE(); i++)
	  if (SurfaceElement(i).GetNP() != 3)
	    return false;
	return true;
      }

    for (int i = 1; i <= GetNSE(); i++)
      if (SurfaceElement(i).GetIndex() == faceindex &&
          SurfaceElement(i).GetNP() != 3)
        return false;
    return true;
  }

  void Mesh ::SetNBCNames ( int nbcn )
  {
    if ( bcnames.Size() )
      for ( int i = 0; i < bcnames.Size(); i++)
        if ( bcnames[i] ) delete bcnames[i];
    bcnames.SetSize(nbcn);
    bcnames = 0;
  }

  void Mesh ::SetBCName ( int bcnr, const string & abcname )
  {
    if ( bcnames[bcnr] ) delete bcnames[bcnr];
    if ( abcname != "default" )
      bcnames[bcnr] = new string ( abcname );
    else
      bcnames[bcnr] = 0;
  }

  const string & Mesh ::GetBCName ( int bcnr ) const
  {
    static string defaultstring = "default";

    if ( !bcnames.Size() )
      return defaultstring;
    if ( bcnames[bcnr] )
      return *bcnames[bcnr];
    else
      return defaultstring;
  }
}
