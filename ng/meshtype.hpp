#pragma once
/**************************************************************************/
/* File:   meshtype.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

namespace netgen
{
  /*
    Classes for NETGEN
  */

  enum ELEMENT_TYPE { 
    SEGMENT = 1, SEGMENT3 = 2,
    TRIG = 10, QUAD=11, TRIG6 = 12, QUAD6 = 13, QUAD8 = 14,
    TET = 20, TET10 = 21, 
    PYRAMID = 22, PRISM = 23, PRISM12 = 24,
    HEX = 25
  };

#define ELEMENT2D_MAXPOINTS 8

  enum POINTTYPE { FIXEDPOINT = 1, EDGEPOINT = 2, SURFACEPOINT = 3, INNERPOINT = 4 };
  enum ELEMENTTYPE { FREEELEMENT, FIXEDELEMENT };
  enum OPTIMIZEGOAL { OPT_QUALITY, OPT_CONFORM, OPT_REST, OPT_WORSTCASE, OPT_LEGAL };

  class PointGeomInfo
  {
  public:
    int trignum;   // for STL Meshing
    double u, v;   // for OCC Meshing

    PointGeomInfo () 
      : trignum(-1), u(0), v(0) { ; }
  };

  inline ostream & operator<< (ostream & ost, const PointGeomInfo & gi)
  {
    return (ost << gi.trignum << " " << gi.u << " " << gi.v);
  }

  inline istream & operator>> (istream & ist, PointGeomInfo & gi)
  {
    return (ist >> gi.trignum >> gi.u >> gi.v);
  }

#define MULTIPOINTGEOMINFO_MAX 100
  class MultiPointGeomInfo
  {
    int cnt;
    PointGeomInfo mgi[MULTIPOINTGEOMINFO_MAX];
  public:
    MultiPointGeomInfo () { cnt = 0; }
    int AddPointGeomInfo (const PointGeomInfo & gi);
    void Init () { cnt = 0; }
    void DeleteAll () { cnt = 0; }

    int GetNPGI () const { return cnt; }
    const PointGeomInfo & GetPGI (int i) const { return mgi[i-1]; }
  };

  class EdgePointGeomInfo
  {
  public:
    int edgenr;
    int body;    // for ACIS
    double dist; // for 2d meshing
    double u, v; // for OCC Meshing

  public:
    EdgePointGeomInfo ()
      : edgenr(0), body(0), dist(0.0), u(0.0), v(0.0) { ; }

    EdgePointGeomInfo & operator= (const EdgePointGeomInfo & gi2)
    {
      edgenr = gi2.edgenr;  
      body = gi2.body;
      dist = gi2.dist;
      u = gi2.u; v = gi2.v;
      return *this;
    }
  };

  inline ostream & operator<< (ostream & ost, const EdgePointGeomInfo & gi)
  {
    ost << "epgi: edgnr=" << gi.edgenr << ", dist=" << gi.dist;
    return ost;
  }

  class PointIndex
  {
    int i;
  public:
    PointIndex () { ; }
    PointIndex (int ai) : i(ai) { ; }
    PointIndex & operator= (const PointIndex &ai) { i = ai.i; return *this; }
    // PointIndex & operator= (int ai) { i = ai; return *this; }
    operator int () const { return i; }
    // int GetInt () const { return i; }
    // PointIndex operator+ (int i2) { return PointIndex (i+i2); }
    // PointIndex operator++ (int) { int hi = i; i++; return PointIndex(hi); }
    // PointIndex operator-- (int) { int hi = i; i--; return PointIndex(hi); }
    PointIndex operator++ (int) { PointIndex hi(*this); i++; return hi; }
    PointIndex operator-- (int) { PointIndex hi(*this); i--; return hi; }
    PointIndex operator++ () { i++; return *this; }
    PointIndex operator-- () { i--; return *this; }

    enum { BASE = 1 };  
  };

  inline istream & operator>> (istream & ist, PointIndex & pi)
  {
    int i; ist >> i; pi = i; return ist;
  }

  inline ostream & operator<< (ostream & ost, const PointIndex & pi)
  {
    return (ost << int(pi));
  }

  class ElementIndex
  {
    int i;
  public:
    ElementIndex () { ; }
    ElementIndex (int ai) : i(ai) { ; }
    ElementIndex & operator= (const ElementIndex & ai) { i = ai.i; return *this; }
    ElementIndex & operator= (int ai) { i = ai; return *this; }
    operator int () const { return i; }
    ElementIndex & operator++ (int) { i++; return *this; }
    ElementIndex & operator-- (int) { i--; return *this; }
  };

  inline istream & operator>> (istream & ist, ElementIndex & pi)
  {
    int i; ist >> i; pi = i; return ist;
  }

  inline ostream & operator<< (ostream & ost, const ElementIndex & pi)
  {
    return (ost << int(pi));
  }

  class SurfaceElementIndex
  {
    int i;
  public:
    SurfaceElementIndex () { ; }
    SurfaceElementIndex (int ai) : i(ai) { ; }
    SurfaceElementIndex & operator= (const SurfaceElementIndex & ai) 
    { i = ai.i; return *this; }
    SurfaceElementIndex & operator= (int ai) { i = ai; return *this; }
    operator int () const { return i; }
    SurfaceElementIndex & operator++ (int) { i++; return *this; }
    SurfaceElementIndex & operator+= (int inc) { i+=inc; return *this; }
    SurfaceElementIndex & operator-- (int) { i--; return *this; }
  };

  inline istream & operator>> (istream & ist, SurfaceElementIndex & pi)
  {
    int i; ist >> i; pi = i; return ist;
  }

  inline ostream & operator<< (ostream & ost, const SurfaceElementIndex & pi)
  {
    return (ost << int(pi));
  }

  class SegmentIndex
  {
    int i;
  public:
    SegmentIndex () { ; }
    SegmentIndex (int ai) : i(ai) { ; }
    SegmentIndex & operator= (const SegmentIndex & ai) 
    { i = ai.i; return *this; }
    SegmentIndex & operator= (int ai) { i = ai; return *this; }
    operator int () const { return i; }
    SegmentIndex & operator++ (int) { i++; return *this; }
    SegmentIndex & operator-- (int) { i--; return *this; }
  };

  inline istream & operator>> (istream & ist, SegmentIndex & pi)
  {
    int i; ist >> i; pi = i; return ist;
  }

  inline ostream & operator<< (ostream & ost, const SegmentIndex & pi)
  {
    return (ost << int(pi));
  }

  /**
     Point in the mesh.
     Contains layer (a new feature in 4.3 for overlapping meshes.
  */
  class MeshPoint : public Point<3>
  {
    int layer;
    double singular; // singular factor for hp-refinement
    POINTTYPE type;

  public:
    MeshPoint () 
    { 
      ;
    }

    MeshPoint (const Point<3> & ap, int alayer = 1, POINTTYPE apt = INNERPOINT)
      : Point<3> (ap), layer(alayer), singular(0.),type(apt) 
    { 
      ;
    }
  
    void SetPoint (const Point<3> & ap)
    { 
      Point<3>::operator= (ap); 
      layer = 0; 
      singular = 0; 
    }

    int GetLayer() const { return layer; }

    POINTTYPE Type() const { return type; }
    void SetType(POINTTYPE at) { type = at; }
 
    double Singularity() const { return singular; }
    void Singularity(double s) { singular = s; }
    bool IsSingular() const { return (singular != 0.0); }

  };

  inline ostream & operator<<(ostream  & s, const MeshPoint & pt)
  { 
    return (s << Point<3> (pt)); 
  }

  typedef Array<MeshPoint, PointIndex::BASE, PointIndex> T_POINTS;

  /**
     Triangle element for surface mesh generation.
  */
  class Element2d
  { 
    /// point numbers
    PointIndex pnum[ELEMENT2D_MAXPOINTS];
    /// geom info of points
    PointGeomInfo geominfo[ELEMENT2D_MAXPOINTS];

    /// surface nr
    int index:16;
    ///
    ELEMENT_TYPE typ:6;
    /// number of points
    unsigned int np:4;
    bool badel:1;
    bool refflag:1;  // marked for refinement
    bool strongrefflag:1;
    bool deleted:1;  // element is deleted

    // Philippose - 08 August 2010
    // Set a new property for each element, to 
    // control whether it is visible or not
    bool visible:1;  // element visible

    /// order for hp-FEM
    unsigned int orderx:6;
    unsigned int ordery:6;

    /// a linked list for all segments in the same face
    SurfaceElementIndex next;

  public:
    ///
    Element2d ();
    ///
    Element2d (int anp);
    ///
     Element2d (ELEMENT_TYPE type);
    ///
    Element2d (int pi1, int pi2, int pi3);
    ///
    Element2d (int pi1, int pi2, int pi3, int pi4);
    ///
    ELEMENT_TYPE GetType () const { return typ; }
    /// 
    void SetType (ELEMENT_TYPE atyp)
    {
      typ = atyp;
      switch (typ)
	{
	case TRIG: np = 3; break;
	case QUAD: np = 4; break;
	case TRIG6: np = 6; break;
	case QUAD6: np = 6; break;
	case QUAD8: np = 8; break;
	default:
	  cout << "Element2d::SetType, illegal type " << typ << endl;
	}
    }
    ///
    int GetNP() const { return np; }
    ///
    int GetNV() const
    {
      switch (typ)
	{
	case TRIG:
	case TRIG6: return 3;

	case QUAD:
	case QUAD8:
	case QUAD6: return 4;
	default:
	    ;
	}
      return np;
    }

    ///
    PointIndex & operator[] (int i) { return pnum[i]; }
    ///
    const PointIndex & operator[] (int i) const { return pnum[i]; }

    FlatArray<const PointIndex> PNums () const 
    { return FlatArray<const PointIndex> (np, &pnum[0]); }
    
    ///
    PointIndex & PNum (int i) { return pnum[i-1]; }
    ///
    const PointIndex & PNum (int i) const { return pnum[i-1]; }
    ///
    PointIndex & PNumMod (int i) { return pnum[(i-1) % np]; }
    ///
    const PointIndex & PNumMod (int i) const { return pnum[(i-1) % np]; }
    ///

    ///
    PointGeomInfo & GeomInfoPi (int i) { return geominfo[i-1]; }
    ///
    const PointGeomInfo & GeomInfoPi (int i) const { return geominfo[i-1]; }
    ///
    PointGeomInfo & GeomInfoPiMod (int i) { return geominfo[(i-1) % np]; }
    ///
    const PointGeomInfo & GeomInfoPiMod (int i) const { return geominfo[(i-1) % np]; }


    void SetIndex (int si) { index = si; }
    ///
    int GetIndex () const { return index; }

    int GetOrder () const { return orderx; }
    void SetOrder (int aorder) { orderx = ordery = aorder; }


    void GetOrder (int & ox, int & oy) const { ox = orderx, oy =ordery;};
    void GetOrder (int & ox, int & oy, int & oz) const { ox = orderx; oy = ordery; oz=0; }
    void SetOrder (int ox, int oy, int  /* oz */) { orderx = ox; ordery = oy;}
    void SetOrder (int ox, int oy) { orderx = ox; ordery = oy;}

    ///
    void GetBox (const T_POINTS & points, Box3d & box) const;
    /// invert orientation
    inline void Invert ();
    ///
    void Invert2 ();
    /// first point number is smallest
    inline void NormalizeNumbering ();
    ///
    void NormalizeNumbering2 ();

    bool BadElement() const { return badel; }

    // friend ostream & operator<<(ostream  & s, const Element2d & el);
    friend class Mesh;

    /// get number of 'integration points'
    int GetNIP () const;
    void GetIntegrationPoint (int ip, Point2d & p, double & weight) const;

    void GetTransformation (int ip, const Array<Point2d> & points,
			    class DenseMatrix & trans) const;
    void GetTransformation (int ip, class DenseMatrix & pmat,
			    class DenseMatrix & trans) const;

    void GetShape (const Point2d & p, class Vector & shape) const;
    void GetShapeNew (const Point<2> & p, class FlatVector & shape) const;
    /// matrix 2 * np
    void GetDShape (const Point2d & p, class DenseMatrix & dshape) const;
    void GetDShapeNew (const Point<2> & p, class MatrixFixWidth<2> & dshape) const;
    /// matrix 2 * np
    void GetPointMatrix (const Array<Point2d> & points,
			 class DenseMatrix & pmat) const; 

    void ComputeIntegrationPointData () const;
  

    double CalcJacobianBadness (const Array<Point2d> & points) const;
    double CalcJacobianBadness (const T_POINTS & points, 
				const Vec<3> & n) const;
    double CalcJacobianBadnessDirDeriv (const Array<Point2d> & points,
					int pi, Vec2d & dir, double & dd) const;

    void Delete () { deleted = 1; pnum[0] = pnum[1] = pnum[2] = pnum[3] = PointIndex::BASE-1; }
    bool IsDeleted () const 
    {
      return deleted; 
    }

    // Philippose - 08 August 2010
    // Access functions for the new property: visible
    void Visible(bool vis = 1) 
    { visible = vis; }
    bool IsVisible () const 
    { return visible; }
   
    void SetRefinementFlag (bool rflag = 1) 
    { refflag = rflag; }
    bool TestRefinementFlag () const
    { return refflag; }

    void SetStrongRefinementFlag (bool rflag = 1) 
    { strongrefflag = rflag; }
    bool TestStrongRefinementFlag () const
    { return strongrefflag; }

  
    SurfaceElementIndex NextElement() { return next; }

    bool operator==(const Element2d & el2) const;

    int HasFace(const Element2d& el) const;
    ///
    int meshdocval;
    ///
    int hp_elnr;

  };

  ostream & operator<<(ostream  & s, const Element2d & el);

  class IntegrationPointData
  {
  public:
    Point<3> p;
    double weight;
    Vector shape;
    DenseMatrix dshape;
  };
  
  /**
     Edge segment.
  */
  class Segment
  {
  public:
    ///
     Segment();
     Segment (const Segment& other);

    ~Segment()
    { ; }

    PointIndex pnums[3];  // p1, p2, pmid

    int edgenr;
    ///
    double singedge_left;
    double singedge_right;

    /// 0.. not first segment of segs, 1..first of class, 2..first of class, inverse
    unsigned int seginfo:2;

    /// surface decoding index
    int si;          
    /// domain number inner side
    int domin;
    /// domain number outer side
    int domout;  
    /// top-level object number of surface
    int tlosurf;
    ///
    PointGeomInfo geominfo[2];

    /// surfaces describing edge
    int surfnr1, surfnr2;
    ///
    EdgePointGeomInfo epgeominfo[2];
    ///
    // int pmid; // for second order
    ///
    int meshdocval;

  private:
    string* bcname;

  public:

    Segment& operator=(const Segment & other);

    int hp_elnr;

    void SetBCName ( string * abcname )
    {
      bcname = abcname;
    }

    string * BCNamePtr () 
    { return bcname; }

    const string * BCNamePtr () const 
    { return bcname; }

    const string & GetBCName () const
    {
      static string defaultstring = "default";
      if (! bcname ) return defaultstring;
      return *bcname;
    }

    int GetNP() const
    {
      return (pnums[2] < 0) ? 2 : 3;
    }

    ELEMENT_TYPE GetType() const
    {
      return (pnums[2] < 0) ? SEGMENT : SEGMENT3;
    }
  
    PointIndex & operator[] (int i) { return pnums[i]; }
    const PointIndex & operator[] (int i) const { return pnums[i]; }

    int GetPartition () const { return 0; }

  };

  ostream & operator<<(ostream  & s, const Segment & seg);

  ///
  class FaceDescriptor
  {
    /// which surface, 0 if not available
    int surfnr;
    /// domain nr inside
    int domin;
    /// domain nr outside
    int domout;
    /// top level object number of surface
    int tlosurf;
    /// boundary condition property
    int bcprop;
    // Philippose - 06/07/2009
    // Add capability to store surface colours along with 
    // other face data
    /// surface colour (Default: R=0.0 ; G=1.0 ; B=0.0)
    Vec3d surfcolour;

    ///
    string * bcname;
    /// root of linked list 
    SurfaceElementIndex firstelement;
  
    double domin_singular;
    double domout_singular;

  public:
     FaceDescriptor();
     FaceDescriptor(int surfnri, int domini, int domouti, int tlosurfi);
     FaceDescriptor(const Segment & seg);
     FaceDescriptor(const FaceDescriptor& other);
     ~FaceDescriptor()  { ; }

     int SegmentFits (const Segment & seg);

    int SurfNr () const { return surfnr; }
    int DomainIn () const { return domin; }
    int DomainOut () const { return domout; }
    int TLOSurface () const { return tlosurf; }
    int BCProperty () const { return bcprop; }


    double DomainInSingular() const { return domin_singular; }
    double DomainOutSingular() const { return domout_singular; }

    // Philippose - 06/07/2009
    // Get Surface colour
    Vec3d SurfColour () const { return surfcolour; }
    const string & GetBCName () const;
    void SetSurfNr (int sn) { surfnr = sn; }
    void SetDomainIn (int di) { domin = di; }
    void SetDomainOut (int dom) { domout = dom; }
    void SetBCProperty (int bc) { bcprop = bc; }
    void SetBCName (string * bcn) { bcname = bcn; }
    // Philippose - 06/07/2009
    // Set the surface colour
    void SetSurfColour (Vec3d colour) { surfcolour = colour; }

    void SetDomainInSingular (double v) { domin_singular = v; }
    void SetDomainOutSingular (double v) { domout_singular = v; }

    SurfaceElementIndex FirstElement() { return firstelement; }
    friend class Mesh;
  };

  ostream & operator<< (ostream  & s, const FaceDescriptor & fd);
 
  class EdgeDescriptor
  {
    int tlosurf;
    int surfnr[2];
  public:
    EdgeDescriptor ()
      : tlosurf(-1)
    { surfnr[0] = surfnr[1] = -1; }

    int SurfNr (int i) const { return surfnr[i]; }
    void SetSurfNr (int i, int nr) { surfnr[i] = nr; }

    int TLOSurface() const { return tlosurf; }
    void SetTLOSurface (int nr) { tlosurf = nr; }
  };

  class MeshingParameters
  {
  public:
    bool enableOutput;
    /// use local h ?
    int uselocalh;
    /// grading for local h
    double grading;
    /// maximal mesh size
    double maxh;
    /// minimal mesh size
    double minh;
    /// check overlapping surfaces (debug)
    int checkoverlap;
    /// check chart boundary (sometimes too restrictive)
    int checkchartboundary;
    /// safty factor for curvatures (elemetns per radius)
    double curvaturesafety;
    /// minimal number of segments per edge
    double segmentsperedge;
    /// weight of element size w.r.t element shape
    double elsizeweight;
    /// init with default values

    /// from mp3:
    /// give up quality class, 2d meshing
    int giveuptol2d;
  
    /// high order element curvature
    int elementorder;
    ///
    int inverttets;
    ///
    int inverttrigs;
    ///
    int autozrefine;

    MeshingParameters ();
  };

  inline void Element2d :: Invert()
  {
    if (typ == TRIG)
      Swap (PNum(2), PNum(3));
    else
      Invert2();
  }

  inline void Element2d :: NormalizeNumbering ()
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
      NormalizeNumbering2();
  }

  /**
     Identification of periodic surfaces, close surfaces, etc. 
  */
  class Identifications
  {
  public:
    enum ID_TYPE { UNDEFINED = 1, PERIODIC = 2, CLOSESURFACES = 3, CLOSEEDGES = 4};
  

  private:
    class Mesh & mesh;

    /// identify points (thin layers, periodic b.c.)  
    INDEX_2_HASHTABLE<int> * identifiedpoints;
  
    /// the same, with info about the id-nr
    INDEX_3_HASHTABLE<int> * identifiedpoints_nr;

    /// sorted by identification nr
    TABLE<INDEX_2> idpoints_table;

    Array<ID_TYPE> type;

    /// number of identifications (or, actually used identifications ?)
    int maxidentnr;

  public:
    ///
     Identifications (class Mesh & amesh);
    ///
     ~Identifications ();

     void Delete ();

    /*
      Identify points pi1 and pi2, due to
      identification nr identnr
    */
     void Add (PointIndex pi1, PointIndex pi2, int identnr);

    int Get (PointIndex pi1, PointIndex pi2) const;
    int GetSymmetric (PointIndex pi1, PointIndex pi2) const;

    bool Get (PointIndex pi1, PointIndex pi2, int identnr) const;
    bool GetSymmetric (PointIndex pi1, PointIndex pi2, int identnr) const;

    ///
    INDEX_2_HASHTABLE<int> & GetIdentifiedPoints () 
    { 
      return *identifiedpoints; 
    }

    bool Used (PointIndex pi1, PointIndex pi2)
    {
      return identifiedpoints->Used (INDEX_2 (pi1, pi2));
    }

    bool UsedSymmetric (PointIndex pi1, PointIndex pi2)
    {
      return 
	identifiedpoints->Used (INDEX_2 (pi1, pi2)) ||
	identifiedpoints->Used (INDEX_2 (pi2, pi1));
    }

    ///
    void GetMap (int identnr, Array<int,PointIndex::BASE> & identmap, bool symmetric = false) const;
    ///
    ID_TYPE GetType(int identnr) const
    {
      if(identnr <= type.Size())
	return type[identnr-1];
      else
	return UNDEFINED;
    }
    void SetType(int identnr, ID_TYPE t)
    {
      while(type.Size() < identnr)
	type.Append(UNDEFINED);
      type[identnr-1] = t;
    }
    
    void GetPairs (int identnr, Array<INDEX_2> & identpairs) const;
    ///
    int GetMaxNr () const { return maxidentnr; }  

    /// remove secondorder
    void SetMaxPointNr (int maxpnum);
  };
}