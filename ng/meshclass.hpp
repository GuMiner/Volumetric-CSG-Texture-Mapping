#pragma once
/**************************************************************************/
/* File:   meshclass.hpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   20. Nov. 99                                                    */
/**************************************************************************/

/*
  The mesh class
*/
namespace netgen
{
  enum resthtype { RESTRICTH_FACE, RESTRICTH_EDGE, 
		   RESTRICTH_SURFACEELEMENT, RESTRICTH_POINT, RESTRICTH_SEGMENT };

  /// 2d/3d mesh
  class Mesh
  {
  public:
    typedef ::netgen::T_POINTS T_POINTS;
    typedef Array<Element2d> T_SURFELEMENTS;

  private:
    /// point coordinates
    T_POINTS points;

    /// line-segments at edges
    Array<Segment> segments;
    /// surface elements, 2d-inner elements
    T_SURFELEMENTS surfelements;
    /// points will be fixed forever
    Array<PointIndex> lockedpoints;

    /// surface indices at boundary nodes
    TABLE<int,PointIndex::BASE> surfacesonnode;
    /// boundary edges  (1..normal bedge, 2..segment)
    INDEX_2_CLOSED_HASHTABLE<int> * boundaryedges;
    ///
    INDEX_2_CLOSED_HASHTABLE<int> * segmentht;
    ///
    INDEX_3_CLOSED_HASHTABLE<int> * surfelementht;

    /// faces of rest-solid
    Array<Element2d> openelements;
    /// open segmenets for surface meshing  
    Array<Segment> opensegments;

    /**
       Representation of local mesh-size h
    */
    LocalH * lochfunc;
    ///
    double hglob;
    ///
    double hmin;
    ///
    Array<double> maxhdomain;
  
    /**
       the face-index of the surface element maps into
       this table.
    */
    Array<FaceDescriptor> facedecoding;

    /**
       the edge-index of the line element maps into
       this table.
    */
    Array<EdgeDescriptor> edgedecoding;

    /// labels for boundary conditions
    Array<string*> bcnames;

    /// Periodic surface, close surface, etc. identifications
    Identifications * ident;

    /// number of vertices (if < 0, use np)
    int numvertices;
  

  private:
    void BuildBoundaryEdges(void);

  public:

    // store coarse mesh before hp-refinement
    Mesh * coarsemesh;

    ///
     Mesh();
    ///
     ~Mesh();

    Mesh & operator= (const Mesh & mesh2);
  
    ///
     void DeleteMesh();
  
    ///
    void ClearSurfaceElements();

    ///
    void ClearSegments()
    { 
      segments.SetSize(0); 
    }

    void SetAllocSize(int nnodes, int nsegs, int nsel, int nel);

     PointIndex AddPoint (const Point3d & p, int layer = 1);
     PointIndex AddPoint (const Point3d & p, int layer, POINTTYPE type);

    int GetNP () const { return points.Size(); }

    MeshPoint & Point(int i) { return points.Elem(i); }
    MeshPoint & Point(PointIndex pi) { return points[pi]; }
    const MeshPoint & Point(int i) const { return points.Get(i); }
    const MeshPoint & Point(PointIndex pi) const { return points[pi]; }

    const MeshPoint & operator[] (PointIndex pi) const { return points[pi]; }
    MeshPoint & operator[] (PointIndex pi) { return points[pi]; }

    const T_POINTS & Points() const { return points; }
    T_POINTS & Points() { return points; }


     SegmentIndex AddSegment (const Segment & s);

    int GetNSeg () const { return segments.Size(); }
    Segment & LineSegment(int i) { return segments.Elem(i); }
    const Segment & LineSegment(int i) const { return segments.Get(i); }

    Segment & LineSegment(SegmentIndex si) { return segments[si]; }
    const Segment & LineSegment(SegmentIndex si) const { return segments[si]; }
    const Segment & operator[] (SegmentIndex si) const { return segments[si]; }
    Segment & operator[] (SegmentIndex si) { return segments[si]; }

     SurfaceElementIndex AddSurfaceElement (const Element2d & el);
    void DeleteSurfaceElement (int eli)
    { 
      surfelements.Elem(eli).Delete();
      surfelements.Elem(eli).PNum(1) = -1; 
      surfelements.Elem(eli).PNum(2) = -1; 
      surfelements.Elem(eli).PNum(3) = -1; 
    }

    void DeleteSurfaceElement (SurfaceElementIndex eli)
    {
      DeleteSurfaceElement (int(eli)+1);
    }

    int GetNSE () const { return surfelements.Size(); }
    Element2d & SurfaceElement(int i) { return surfelements.Elem(i); }
    const Element2d & SurfaceElement(int i) const { return surfelements.Get(i); }
    Element2d & SurfaceElement(SurfaceElementIndex i) { return surfelements[i]; }
    const Element2d & SurfaceElement(SurfaceElementIndex i) const { return surfelements[i]; }

    const Element2d & operator[] (SurfaceElementIndex ei) const
    { return surfelements[ei]; }
    Element2d & operator[] (SurfaceElementIndex ei)
    { return surfelements[ei]; }

     void RebuildSurfaceElementLists ();
     void GetSurfaceElementsOfFace (int facenr, Array<SurfaceElementIndex> & sei) const;

    /// 
     void AddLockedPoint (PointIndex pi);
    ///
    void ClearLockedPoints ();

    const Array<PointIndex> & LockedPoints() const
    { return lockedpoints; }

    /// sets internal tables
    void CalcSurfacesOfNode ();

    /// additional (temporarily) fix points 
    void FixPoints (const BitArray & fixpoints);

    int GetNOpenSegments () { return opensegments.Size(); }
    const Segment & GetOpenSegment (int nr) { return opensegments.Get(nr); }

    /**
       finds average h of surface surfnr if surfnr > 0,
       else of all surfaces.
    */
     double AverageH (int surfnr = 0) const;
    /// Calculates localh 
     void CalcLocalH (double grading);
    ///
     void SetLocalH (const Point3d & pmin, const Point3d & pmax, double grading);
    ///
     void RestrictLocalH (const Point3d & p, double hloc);
    ///
     void RestrictLocalHLine (const Point3d & p1, const Point3d & p2, 
			     double hloc);
    /// number of elements per radius
    ///
     void RestrictLocalH (resthtype rht, int nr, double loch);
    ///
     void SetGlobalH (double h);
    ///
    void SetMinimalH (double h);
    ///
    double MaxHDomain (int dom) const;
    ///
    void SetMaxHDomain (const Array<double> & mhd);
    ///
    double GetH (const Point3d & p) const;
    ///
    double GetMinH (const Point3d & pmin, const Point3d & pmax);
    ///
    LocalH & LocalHFunction () { return * lochfunc; }

    /// Find bounding box
     void GetBox (Point3d & pmin, Point3d & pmax, int dom = -1) const;

    /// Find bounding box of points of typ ptyp or less
     void GetBox (Point3d & pmin, Point3d & pmax, POINTTYPE ptyp ) const;

    ///
    int GetNOpenElements() const
    { return openelements.Size(); }
    ///
    const Element2d & OpenElement(int i) const
    { return openelements.Get(i); }

    /// 
    void SplitSeparatedFaces ();

    bool BoundaryEdge (PointIndex pi1, PointIndex pi2) const
    {
      if(!boundaryedges)
	const_cast<Mesh *>(this)->BuildBoundaryEdges();

      INDEX_2 i2 (pi1, pi2);
      i2.Sort();
      return boundaryedges->Used (i2);
    }

    bool IsSegment (PointIndex pi1, PointIndex pi2) const
    {
      INDEX_2 i2 (pi1, pi2);
      i2.Sort();
      return segmentht->Used (i2);
    }

    SegmentIndex SegmentNr (PointIndex pi1, PointIndex pi2) const
    {
      INDEX_2 i2 (pi1, pi2);
      i2.Sort();
      return segmentht->Get (i2);
    }

    /**
       Remove unused points. etc.
    */
     void Compress ();

    ///
    void ImproveMesh (const MeshingParameters & mp, OPTIMIZEGOAL goal = OPT_QUALITY);

    ///
    bool LegalTrig (const Element2d & el) const;
    
    /// orient surface mesh, for one sub-domain only
    void SurfaceMeshOrientation ();

    ///
    int AddFaceDescriptor(const FaceDescriptor& fd)
    { return facedecoding.Append(fd); }

    int AddEdgeDescriptor(const EdgeDescriptor & fd)
    { return edgedecoding.Append(fd) - 1; }

    ///
     void SetMaterial (int domnr, const char * mat);
    ///
    const char * GetMaterial (int domnr) const;
    
     void SetNBCNames ( int nbcn );

     void SetBCName ( int bcnr, const string & abcname );

    const string & GetBCName ( int bcnr ) const;

    string * GetBCNamePtr ( int bcnr )
    { return bcnames[bcnr]; }

    ///
    void ClearFaceDescriptors()
    { facedecoding.SetSize(0); }

    ///
    int GetNFD () const
    { return facedecoding.Size(); }

    const FaceDescriptor & GetFaceDescriptor (int i) const
    { return facedecoding.Get(i); }

    const EdgeDescriptor & GetEdgeDescriptor (int i) const
    { return edgedecoding[i]; }

    ///
    FaceDescriptor & GetFaceDescriptor (int i) 
    { return facedecoding.Elem(i); }

    /// return periodic, close surface etc. identifications
    Identifications & GetIdentifications () { return *ident; }
    /// return periodic, close surface etc. identifications
    const Identifications & GetIdentifications () const { return *ident; }

    /// find number of vertices
    void ComputeNVertices ();
    /// number of vertices (no edge-midpoints)
    int GetNV () const;
    
    bool PureTrigMesh (int faceindex = 0) const;

    class CSurfaceArea
    {
      const Mesh & mesh;
      bool valid;
      double area;
    public:
      CSurfaceArea (const Mesh & amesh) 
	: mesh(amesh), valid(false) { ; }

      void Add (const Element2d & sel)
      {
	if (sel.GetNP() == 3)
	  area += Cross ( mesh[sel[1]]-mesh[sel[0]],
			  mesh[sel[2]]-mesh[sel[0]] ).Length() / 2;
	else
	  area += Cross (Vec3d (mesh.Point (sel.PNum(1)),
				mesh.Point (sel.PNum(3))),
			 Vec3d (mesh.Point (sel.PNum(1)),
				mesh.Point (sel.PNum(4)))).Length() / 2;;
      }
      void ReCalc ()
      {
	area = 0;
	for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
	  Add (mesh[sei]);
	valid = true;
      }

      operator double () const { return area; }
      bool Valid() const { return valid; }
    };

    CSurfaceArea surfarea;
    CSurfaceArea & SurfaceArea() { return surfarea; }
    const CSurfaceArea & SurfaceArea() const { return surfarea; }

    ///
    friend void OptimizeRestart (Mesh & mesh3d);
    ///

    enum GEOM_TYPE { NO_GEOM = 0, GEOM_2D = 1, GEOM_CSG = 10, GEOM_STL = 11, GEOM_OCC = 12, GEOM_ACIS = 13 };
    GEOM_TYPE geomtype;
  };
}