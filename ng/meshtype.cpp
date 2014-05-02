#include "stdafx.h" 
#include "mystdlib.h"

#include "meshing.hpp"  

namespace netgen
{
  int MultiPointGeomInfo :: 
  AddPointGeomInfo (const PointGeomInfo & gi)
  {
    for (int k = 0; k < cnt; k++)
      if (mgi[k].trignum == gi.trignum)
	return 0;
  
    if (cnt < MULTIPOINTGEOMINFO_MAX)
      {
	mgi[cnt] = gi;
	cnt++;
	return 0;
      }

    throw NgException ("Please report error: MPGI Size too small\n");
  }

  Segment :: Segment() 
  {
    pnums[0] = -1;
    pnums[1] = -1; 
    edgenr = -1;

    singedge_left = 0.;
    singedge_right = 0.;
    seginfo = 0;

    si = -1;

    domin = -1;
    domout = -1;
    tlosurf = -1;

    surfnr1 = -1;
    surfnr2 = -1;
    pnums[2] = -1;
    meshdocval = 0;
    bcname = 0;
  }    

  Segment::Segment (const Segment & other)
    :
    edgenr(other.edgenr),
    singedge_left(other.singedge_left),
    singedge_right(other.singedge_right),
    seginfo(other.seginfo),
    si(other.si),
    domin(other.domin),
    domout(other.domout),
    tlosurf(other.tlosurf),
    geominfo(),
    surfnr1(other.surfnr1),
    surfnr2(other.surfnr2),
    epgeominfo(),
    meshdocval(other.meshdocval),
    hp_elnr(other.hp_elnr)
  {
    for (int j = 0; j < 3; j++)
      pnums[j] = other.pnums[j];

    geominfo[0] = other.geominfo[0];
    geominfo[1] = other.geominfo[1];
    epgeominfo[0] = other.epgeominfo[0];
    epgeominfo[1] = other.epgeominfo[1];
    bcname = other.bcname;
  }

  Segment& Segment::operator=(const Segment & other)
  {
    if (&other != this)
      {
	pnums[0] = other[0];
	pnums[1] = other[1];
	edgenr = other.edgenr;
	singedge_left = other.singedge_left;
	singedge_right = other.singedge_right;
	seginfo = other.seginfo;
	si = other.si;
	domin = other.domin;
	domout = other.domout;
	tlosurf = other.tlosurf;
	geominfo[0] = other.geominfo[0];
	geominfo[1] = other.geominfo[1];
	surfnr1 = other.surfnr1;
	surfnr2 = other.surfnr2;
	epgeominfo[0] = other.epgeominfo[0];
	epgeominfo[1] = other.epgeominfo[1];
	pnums[2] = other.pnums[2];
	meshdocval = other.meshdocval;
	hp_elnr = other.hp_elnr;
	bcname = other.bcname;
      }
    
    return *this;
  }

  ostream & operator<<(ostream  & s, const Segment & seg)
  {
    s << seg[0] << "(gi=" << seg.geominfo[0].trignum << ") - "
      << seg[1] << "(gi=" << seg.geominfo[1].trignum << ")"
      << " domin = " << seg.domin << ", domout = " << seg.domout 
      << " si = " << seg.si << ", edgenr = " << seg.edgenr;
    return s;
  }

  Element2d :: Element2d ()
  { 
    for (int i = 0; i < ELEMENT2D_MAXPOINTS; i++)
      {
	pnum[i] = 0;
	geominfo[i].trignum = 0;
      }
    np = 3;
    index = 0;
    badel = 0;
    deleted = 0;
    visible = 1;
    typ = TRIG;
    orderx = ordery = 1;
    refflag = 1;
    strongrefflag = false;
  } 

  Element2d :: Element2d (int anp)
  { 
    for (int i = 0; i < ELEMENT2D_MAXPOINTS; i++)
      {
	pnum[i] = 0;
	geominfo[i].trignum = 0;
      }
    np = anp;
    index = 0;
    badel = 0;
    deleted = 0;
    visible = 1;
    switch (np)
      {
      case 3: typ = TRIG; break;
      case 4: typ = QUAD; break;
      case 6: typ = TRIG6; break;
      case 8: typ = QUAD8; break;
      }
    orderx = ordery = 1;
    refflag = 1;
    strongrefflag = false;
  } 

  Element2d :: Element2d (ELEMENT_TYPE atyp)
  { 
    for (int i = 0; i < ELEMENT2D_MAXPOINTS; i++)
      {
	pnum[i] = 0;
	geominfo[i].trignum = 0;
      }

    SetType (atyp);

    index = 0;
    badel = 0;
    deleted = 0;
    visible = 1;
    orderx = ordery = 1;
    refflag = 1;
    strongrefflag = false;
  } 

  Element2d :: Element2d (int pi1, int pi2, int pi3)
  {
    pnum[0] = pi1;
    pnum[1] = pi2;
    pnum[2] = pi3;
    np = 3;
    typ = TRIG;
    pnum[3] = 0;
    pnum[4] = 0;
    pnum[5] = 0;
  
    for (int i = 0; i < ELEMENT2D_MAXPOINTS; i++)
      geominfo[i].trignum = 0;
    index = 0;
    badel = 0;
    refflag = 1;
    strongrefflag = false;
    deleted = 0;
    visible = 1;
    orderx = ordery = 1;
  }

  Element2d :: Element2d (int pi1, int pi2, int pi3, int pi4)
  {
    pnum[0] = pi1;
    pnum[1] = pi2;
    pnum[2] = pi3;
    pnum[3] = pi4;
    np = 4;
    typ = QUAD;

    pnum[4] = 0;
    pnum[5] = 0;
  
    for (int i = 0; i < ELEMENT2D_MAXPOINTS; i++)
      geominfo[i].trignum = 0;
    index = 0;
    badel = 0;
    refflag = 1;
    strongrefflag = false;
    deleted = 0;
    visible = 1;
    orderx = ordery = 1;
  }

  void Element2d :: GetBox (const T_POINTS & points, Box3d & box) const
  {
    box.SetPoint (points.Get(pnum[0]));
    for (unsigned i = 1; i < np; i++)
      box.AddPoint (points.Get(pnum[i]));
  }

  bool Element2d :: operator==(const Element2d & el2) const
  {
    bool retval = (el2.GetNP() == np);
    for(int i= 0; retval && i<(int)np; i++)
      retval = (el2[i] == (*this)[i]);

    return retval;
  }

  void Element2d :: Invert2()
  {
    switch (typ)
      {
      case TRIG:
        {
          Swap (pnum[1], pnum[2]);
          break;
        }
      case TRIG6:
        {
          Swap (pnum[1], pnum[2]);
          Swap (pnum[4], pnum[5]);
          break;
        }
      case QUAD:
        {
          Swap (pnum[0], pnum[3]);
          Swap (pnum[1], pnum[2]);
          break;
        }
      default:
        {
          cerr << "Element2d::Invert2, illegal element type " << int(typ) << endl;
        }
      }
  }

  int Element2d::HasFace(const Element2d & el) const
  {
    //nur für tets!!! hannes
    for (int i = 1; i <= 3; i++)
      {
        if (PNumMod(i)   == el[0] && 
            PNumMod(i+1) == el[1] && 
            PNumMod(i+2) == el[2])
          {
            return 1;
          }
      }
    return 0;
  }

  void Element2d :: NormalizeNumbering2 ()
  {
    if (GetNP() == 3)
      {
        if (PNum(1) < PNum(2) && PNum(1) < PNum(3))
          return;
        else
          {
            if (PNum(2) < PNum(3))
              {
                PointIndex pi1 = PNum(2);
                PNum(2) = PNum(3);
                PNum(3) = PNum(1);
                PNum(1) = pi1;
              }
            else
              {
                PointIndex pi1 = PNum(3);
                PNum(3) = PNum(2);
                PNum(2) = PNum(1);
                PNum(1) = pi1;
              }
          }
      }
    else
      {
        int mini = 1;
        for (int i = 2; i <= GetNP(); i++)
          if (PNum(i) < PNum(mini)) mini = i;
      
        Element2d hel = (*this);
        for (int i = 1; i <= GetNP(); i++)
          PNum(i) = hel.PNumMod (i+mini-1);
      }
  }

  Array<IntegrationPointData*> ipdtrig;
  Array<IntegrationPointData*> ipdquad;

  int Element2d :: GetNIP () const
  {
    int nip;
    switch (np)
      {
      case 3: nip = 1; break;
      case 4: nip = 4; break;
      default: nip = 0; break;
      }
    return nip;
  }

  void Element2d :: 
  GetIntegrationPoint (int ip, Point2d & p, double & weight) const
  {
    static double eltriqp[1][3] =
      {
        { 1.0/3.0, 1.0/3.0, 0.5 }
      };

    static double elquadqp[4][3] =
      { 
        { 0, 0, 0.25 },
        { 0, 1, 0.25 },
        { 1, 0, 0.25 },
        { 1, 1, 0.25 }
      };
  
    double * pp = 0;
    switch (typ)
      {
      case TRIG: pp = &eltriqp[0][0]; break;
      case QUAD: pp = &elquadqp[ip-1][0]; break;
      default:
        cout << "Element2d::GetIntegrationPoint, illegal type " << typ << std::endl;
      }

    p.X() = pp[0];
    p.Y() = pp[1];
    weight = pp[2];
  }

  void Element2d :: 
  GetTransformation (int ip, const Array<Point2d> & points,
                     DenseMatrix & trans) const
  {
    int np = GetNP();
    DenseMatrix pmat(2, np), dshape(2, np);
    pmat.SetSize (2, np);
    dshape.SetSize (2, np);

    Point2d p;
    double w;

    GetPointMatrix (points, pmat);
    GetIntegrationPoint (ip, p, w);
    GetDShape (p, dshape);
  
    CalcABt (pmat, dshape, trans);
  }

  void Element2d :: 
  GetTransformation (int ip, class DenseMatrix & pmat,
                     class DenseMatrix & trans) const
  {
    ComputeIntegrationPointData ();
    DenseMatrix * dshapep = NULL;
    switch (typ)
      {
      case TRIG: dshapep = &ipdtrig.Get(ip)->dshape; break;
      case QUAD: dshapep = &ipdquad.Get(ip)->dshape; break;
      default:
        cout << "Element2d::GetTransformation, illegal type " << typ << endl;
      }
  
    CalcABt (pmat, *dshapep, trans);
  }

  void Element2d :: GetShape (const Point2d & p, Vector & shape) const
  {
    if (shape.Size() != GetNP())
      {
        cerr << "Element::GetShape: Length not fitting" << endl;
        return;
      }

    switch (typ)
      {
      case TRIG:
        shape(0) = 1 - p.X() - p.Y();
        shape(1) = p.X();
        shape(2) = p.Y();
        break;
      case QUAD:
        shape(0) = (1-p.X()) * (1-p.Y());
        shape(1) = p.X() * (1-p.Y());
        shape(2) = p.X() * p.Y();
        shape(3) = (1-p.X()) * p.Y();
        break;
      default:
        cout <<"Element2d::GetShape, illegal type " << typ << endl;
      }
  }

  void Element2d :: GetShapeNew (const Point<2> & p, FlatVector & shape) const
  {
    switch (typ)
      {
      case TRIG:
        {
          shape(0) = p(0);
          shape(1) = p(1);
          shape(2) = 1-p(0)-p(1);
          break;
        }

      case QUAD:
        {
          shape(0) = (1-p(0))*(1-p(1));
          shape(1) =    p(0) *(1-p(1));
          shape(2) =    p(0) *   p(1) ;
          shape(3) = (1-p(0))*   p(1) ;
          break;
        }
      }
  }

  void Element2d :: 
  GetDShape (const Point2d & p, DenseMatrix & dshape) const
  {
    switch (typ)
      {
      case TRIG:
        dshape.Elem(1, 1) = -1;
        dshape.Elem(1, 2) = 1;
        dshape.Elem(1, 3) = 0;
        dshape.Elem(2, 1) = -1;
        dshape.Elem(2, 2) = 0;
        dshape.Elem(2, 3) = 1;
        break;
      case QUAD:
        dshape.Elem(1, 1) = -(1-p.Y());
        dshape.Elem(1, 2) = (1-p.Y());
        dshape.Elem(1, 3) = p.Y();
        dshape.Elem(1, 4) = -p.Y();
        dshape.Elem(2, 1) = -(1-p.X());
        dshape.Elem(2, 2) = -p.X();
        dshape.Elem(2, 3) = p.X();
        dshape.Elem(2, 4) = (1-p.X());
        break;

      default:
        cout << "Element2d::GetDShape, illegal type " << typ << endl;
      }
  }

  void Element2d :: 
  GetDShapeNew (const Point<2> & p, MatrixFixWidth<2> & dshape) const
  {
    switch (typ)
      {
      case TRIG:
        {
          dshape = 0;
          dshape(0,0) = 1;
          dshape(1,1) = 1;
          dshape(2,0) = -1;
          dshape(2,1) = -1;
          break;
        }
      case QUAD:
        {
          dshape(0,0) = -(1-p(1));
          dshape(0,1) = -(1-p(0));

          dshape(1,0) =  (1-p(1));
          dshape(1,1) =  -p(0);

          dshape(2,0) = p(1);
          dshape(2,1) = p(0);

          dshape(3,0) = -p(1);
          dshape(3,1) = (1-p(0));
          break;
        }
      }
  }

  void Element2d :: 
  GetPointMatrix (const Array<Point2d> & points,
                  DenseMatrix & pmat) const
  {
    int np = GetNP();
  
    for (int i = 1; i <= np; i++)
      {
        const Point2d & p = points.Get(PNum(i));
        pmat.Elem(1, i) = p.X();
        pmat.Elem(2, i) = p.Y();
      }
  }

  double Element2d :: CalcJacobianBadness (const Array<Point2d> & points) const
  {
    int i, j;
    int nip = GetNIP();
    DenseMatrix trans(2,2);
    DenseMatrix pmat;
  
    pmat.SetSize (2, GetNP());
    GetPointMatrix (points, pmat);

    double err = 0;
    for (i = 1; i <= nip; i++)
      {
        GetTransformation (i, pmat, trans);

        // Frobenius norm
        double frob = 0;
        for (j = 1; j <= 4; j++)
          frob += sqr (trans.Get(j));
        frob = sqrt (frob);
        frob /= 2;

        double det = trans.Det();

        if (det <= 0)
          err += 1e12;
        else
          err += frob * frob / det;
      }

    err /= nip;
    return err;
  }
  
  static const int qip_table[4][4] =
    { { 0, 1, 0, 3 },
      { 0, 1, 1, 2 },
      { 3, 2, 0, 3 },
      { 3, 2, 1, 2 }
    };

  double Element2d :: 
  CalcJacobianBadnessDirDeriv (const Array<Point2d> & points,
                               int pi, Vec2d & dir, double & dd) const
  {
    if (typ == QUAD)
      {
        Mat<2,2> trans, dtrans;
        Mat<2,4> vmat, pmat;
      
        for (int j = 0; j < 4; j++)
          {
            const Point2d & p = points.Get( (*this)[j] );
            pmat(0, j) = p.X();
            pmat(1, j) = p.Y();
          }

        vmat = 0.0;
        vmat(0, pi-1) = dir.X();
        vmat(1, pi-1) = dir.Y();
      
        double err = 0;
        dd = 0;

        for (int i = 0; i < 4; i++)
          {
            int ix1 = qip_table[i][0];
            int ix2 = qip_table[i][1];
            int iy1 = qip_table[i][2];
            int iy2 = qip_table[i][3];
	      
            trans(0,0) = pmat(0, ix2) - pmat(0,ix1);
            trans(1,0) = pmat(1, ix2) - pmat(1,ix1);
            trans(0,1) = pmat(0, iy2) - pmat(0,iy1);
            trans(1,1) = pmat(1, iy2) - pmat(1,iy1);

            double det = trans(0,0)*trans(1,1)-trans(1,0)*trans(0,1);

            if (det <= 0)
              {
                dd = 0;
                return 1e12;
              }
	  
            dtrans(0,0) = vmat(0, ix2) - vmat(0,ix1);
            dtrans(1,0) = vmat(1, ix2) - vmat(1,ix1);
            dtrans(0,1) = vmat(0, iy2) - vmat(0,iy1);
            dtrans(1,1) = vmat(1, iy2) - vmat(1,iy1);


            // Frobenius norm
            double frob = 0;
            for (int j = 0; j < 4; j++) 
              frob += sqr (trans(j));
            frob = sqrt (frob);
	  
            double dfrob = 0;
            for (int j = 0; j < 4; j++)
              dfrob += trans(j) * dtrans(j);
            dfrob = dfrob / frob;
	  
            frob /= 2;      
            dfrob /= 2;
	  
            // ddet = \sum_j det (m_j)   with m_j = trans, except col j = dtrans
            double ddet 
              = dtrans(0,0) * trans(1,1) - trans(0,1) * dtrans(1,0)
              + trans(0,0) * dtrans(1,1) - dtrans(0,1) * trans(1,0);
	  
            err += frob * frob / det;
            dd += (2 * frob * dfrob * det - frob * frob * ddet) / (det * det);
          }
      
        err /= 4;
        dd /= 4;
        return err;
      }

    int nip = GetNIP();
    DenseMatrix trans(2,2), dtrans(2,2);
    DenseMatrix pmat, vmat;
  
    pmat.SetSize (2, GetNP());
    vmat.SetSize (2, GetNP());

    GetPointMatrix (points, pmat);
  
    vmat = 0.0;
    vmat.Elem(1, pi) = dir.X();
    vmat.Elem(2, pi) = dir.Y();

    double err = 0;
    dd = 0;

    for (int i = 1; i <= nip; i++)
      {
        GetTransformation (i, pmat, trans);
        GetTransformation (i, vmat, dtrans);

        // Frobenius norm
        double frob = 0;
        for (int j = 1; j <= 4; j++)
          frob += sqr (trans.Get(j));
        frob = sqrt (frob);
      
        double dfrob = 0;
        for (int j = 1; j <= 4; j++)
          dfrob += trans.Get(j) * dtrans.Get(j);
        dfrob = dfrob / frob;
      
        frob /= 2;      
        dfrob /= 2;
      
        double det = trans(0,0)*trans(1,1)-trans(1,0)*trans(0,1);

        // ddet = \sum_j det (m_j)   with m_j = trans, except col j = dtrans
        double ddet 
          = dtrans(0,0) * trans(1,1) - trans(0,1) * dtrans(1,0)
          + trans(0,0) * dtrans(1,1) - dtrans(0,1) * trans(1,0);

        if (det <= 0)
          err += 1e12;
        else
          {
            err += frob * frob / det;
            dd += (2 * frob * dfrob * det - frob * frob * ddet) / (det * det);
          }
      }

    err /= nip;
    dd /= nip;
    return err;
  }

  double Element2d :: 
  CalcJacobianBadness (const T_POINTS & points, const Vec<3> & n) const
  {
    int i, j;
    int nip = GetNIP();
    DenseMatrix trans(2,2);
    DenseMatrix pmat;
  
    pmat.SetSize (2, GetNP());

    Vec<3> t1, t2;
    t1 = n.GetNormal();
    t2 = Cross (n, t1);

    for (i = 1; i <= GetNP(); i++)
      {
        Point3d p = points.Get(PNum(i));
        pmat.Elem(1, i) = p.X() * t1(0) + p.Y() * t1(1) + p.Z() * t1(2);
        pmat.Elem(2, i) = p.X() * t2(0) + p.Y() * t2(1) + p.Z() * t2(2);
      }

    double err = 0;
    for (i = 1; i <= nip; i++)
      {
        GetTransformation (i, pmat, trans);

        // Frobenius norm
        double frob = 0;
        for (j = 1; j <= 4; j++)
          frob += sqr (trans.Get(j));
        frob = sqrt (frob);
        frob /= 2;

        double det = trans.Det();
        if (det <= 0)
          err += 1e12;
        else
          err += frob * frob / det;
      }

    err /= nip;
    return err;
  }

  void Element2d :: ComputeIntegrationPointData () const
  {
    switch (np)
      {
      case 3: if (ipdtrig.Size()) return; break;
      case 4: if (ipdquad.Size()) return; break;
      }

    for (int i = 1; i <= GetNIP(); i++)
      {
        IntegrationPointData * ipd = new IntegrationPointData;
        Point2d hp;
        GetIntegrationPoint (i, hp, ipd->weight);
        ipd->p(0) = hp.X();
        ipd->p(1) = hp.Y();
        ipd->p(2) = 0;

        ipd->shape.SetSize(GetNP());
        ipd->dshape.SetSize(2, GetNP());

        GetShape (hp, ipd->shape);
        GetDShape (hp, ipd->dshape);

        switch (np)
          {
          case 3: ipdtrig.Append (ipd); break;
          case 4: ipdquad.Append (ipd); break;
          }
      }
  }

  ostream & operator<<(ostream  & s, const Element2d & el)
  {
    s << "np = " << el.GetNP();
    for (int j = 1; j <= el.GetNP(); j++)
      s << " " << el.PNum(j);
    return s;
  }

  FaceDescriptor ::  FaceDescriptor()
  { 
    surfnr = domin = domout  = bcprop = 0; 
    domin_singular = domout_singular = 0.;
    // Philippose - 06/07/2009
    // Initialise surface colour
    surfcolour = Vec3d(0.0,1.0,0.0);
    tlosurf = -1; 
    bcname = 0;
    firstelement = -1;
  }

  FaceDescriptor ::  FaceDescriptor(const FaceDescriptor& other)
    : surfnr(other.surfnr), domin(other.domin), domout(other.domout),
      tlosurf(other.tlosurf), bcprop(other.bcprop), 
      surfcolour(other.surfcolour), bcname(other.bcname),
      domin_singular(other.domin_singular), domout_singular(other.domout_singular)
  { 
    firstelement = -1;
  }

  FaceDescriptor :: 
  FaceDescriptor(int surfnri, int domini, int domouti, int tlosurfi)
  { 
    surfnr = surfnri; 
    domin = domini; 
    domout = domouti;
    // Philippose - 06/07/2009
    // Initialise surface colour
    surfcolour = Vec3d(0.0,1.0,0.0);
    tlosurf = tlosurfi; 
    bcprop = surfnri;
    domin_singular = domout_singular = 0.;
    bcname = 0;
    firstelement = -1;
  }

  FaceDescriptor :: FaceDescriptor(const Segment & seg)
  { 
    surfnr = seg.si; 
    domin = seg.domin+1;
    domout = seg.domout+1;
    // Philippose - 06/07/2009
    // Initialise surface colour
    surfcolour = Vec3d(0.0,1.0,0.0);
    tlosurf = seg.tlosurf+1;
    bcprop = 0;
    domin_singular = domout_singular = 0.;
    bcname = 0;
    firstelement = -1;
  }

  int FaceDescriptor ::  SegmentFits (const Segment & seg)
  {
    return
      surfnr == seg.si &&
      domin == seg.domin+1 &&
      domout == seg.domout+1  &&
      tlosurf == seg.tlosurf+1;
  }

  const string & FaceDescriptor :: GetBCName () const
  {
    static string defaultstring = "default";
    if (bcname) return *bcname;
    return defaultstring;
  }

  ostream & operator<<(ostream  & s, const FaceDescriptor & fd)
  {
    s << "surfnr = " << fd.SurfNr() 
      << ", domin = " << fd.DomainIn()
      << ", domout = " << fd.DomainOut()
      << ", tlosurf = " << fd.TLOSurface()
      << ", bcprop = " << fd.BCProperty()
      << ", domin_sing = " << fd.DomainInSingular()
      << ", domout_sing = " << fd.DomainOutSingular()
      << ", colour = " << fd.SurfColour();
    return s;
  }

  Identifications :: Identifications (Mesh & amesh)
    : mesh(amesh)
  {
    identifiedpoints = new INDEX_2_HASHTABLE<int>(100);
    identifiedpoints_nr = new INDEX_3_HASHTABLE<int>(100);
    maxidentnr = 0;
  }

  Identifications :: ~Identifications ()
  {
    delete identifiedpoints;
    delete identifiedpoints_nr;
  }

  void Identifications :: Delete ()
  {
    delete identifiedpoints;
    identifiedpoints = new INDEX_2_HASHTABLE<int>(100);
    delete identifiedpoints_nr;
    identifiedpoints_nr = new INDEX_3_HASHTABLE<int>(100);
    maxidentnr = 0;
  }

  void Identifications :: Add (PointIndex pi1, PointIndex pi2, int identnr)
  {
    INDEX_2 pair (pi1, pi2);
    identifiedpoints->Set (pair, identnr);

    INDEX_3 tripl (pi1, pi2, identnr);
    identifiedpoints_nr->Set (tripl, 1);

    if (identnr > maxidentnr) maxidentnr = identnr;

    if (identnr+1 > idpoints_table.Size())
      idpoints_table.ChangeSize (identnr+1);
    idpoints_table.Add (identnr, pair);
  }

  int Identifications :: Get (PointIndex pi1, PointIndex pi2) const
  {
    INDEX_2 pair(pi1, pi2);
    if (identifiedpoints->Used (pair))
      return identifiedpoints->Get(pair);
    else
      return 0;
  }

  bool Identifications :: Get (PointIndex pi1, PointIndex pi2, int nr) const
  {
    INDEX_3 tripl(pi1, pi2, nr);
    if (identifiedpoints_nr->Used (tripl))
      return 1;
    else
      return 0;
  }

  int Identifications :: GetSymmetric (PointIndex pi1, PointIndex pi2) const
  {
    INDEX_2 pair(pi1, pi2);
    if (identifiedpoints->Used (pair))
      return identifiedpoints->Get(pair);

    pair = INDEX_2 (pi2, pi1);
    if (identifiedpoints->Used (pair))
      return identifiedpoints->Get(pair);

    return 0;
  }


  void Identifications :: GetMap (int identnr, Array<int,PointIndex::BASE> & identmap, bool symmetric) const
  {
    identmap.SetSize (mesh.GetNP());
    identmap = 0;

    if (identnr)
      for (int i = 0; i < idpoints_table[identnr].Size(); i++)
        {
          INDEX_2 pair = idpoints_table[identnr][i];
          identmap[pair.I1()] = pair.I2();
          if(symmetric)
            identmap[pair.I2()] = pair.I1();
        }

    else
      {

        for (int i = 1; i <= identifiedpoints_nr->GetNBags(); i++)
          for (int j = 1; j <= identifiedpoints_nr->GetBagSize(i); j++)
            {
              INDEX_3 i3;
              int dummy;
              identifiedpoints_nr->GetData (i, j, i3, dummy);
	    
              if (i3.I3() == identnr || !identnr)
                {
                  identmap.Elem(i3.I1()) = i3.I2();
                  if(symmetric)
                    identmap.Elem(i3.I2()) = i3.I1();
                }
            }  
      }

  }

  void Identifications :: GetPairs (int identnr, 
                                    Array<INDEX_2> & identpairs) const
  {
    identpairs.SetSize(0);
  
    if (identnr == 0)
      for (int i = 1; i <= identifiedpoints->GetNBags(); i++)
        for (int j = 1; j <= identifiedpoints->GetBagSize(i); j++)
          {
            INDEX_2 i2;
            int nr;
            identifiedpoints->GetData (i, j, i2, nr);
            identpairs.Append (i2);
          }  
    else
      for (int i = 1; i <= identifiedpoints_nr->GetNBags(); i++)
        for (int j = 1; j <= identifiedpoints_nr->GetBagSize(i); j++)
          {
            INDEX_3 i3;
            int dummy;
            identifiedpoints_nr->GetData (i, j, i3 , dummy);
	  
            if (i3.I3() == identnr)
              identpairs.Append (INDEX_2(i3.I1(), i3.I2()));
          }  
  }

  void Identifications :: SetMaxPointNr (int maxpnum)
  {
    for (int i = 1; i <= identifiedpoints->GetNBags(); i++)
      for (int j = 1; j <= identifiedpoints->GetBagSize(i); j++)
        {
          INDEX_2 i2;
          int nr;
          identifiedpoints->GetData (i, j, i2, nr);
	
          if (i2.I1() > maxpnum || i2.I2() > maxpnum)
            {
              i2.I1() = i2.I2() = -1;
              identifiedpoints->SetData (i, j, i2, -1);	    
            }
        }
  }

  MeshingParameters :: MeshingParameters ()
  {
    enableOutput = false;
    uselocalh = 1;
    grading = 0.3;
    maxh = 1e10;
    minh = 0;
    checkoverlap = 1;
    checkchartboundary = 1;
    curvaturesafety = 2;
    segmentsperedge = 1;
    elsizeweight = 0.2;
    giveuptol2d = 200;
  }
}