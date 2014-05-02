#include "stdafx.h" 
#include "mystdlib.h"
#include "meshing.hpp"

namespace netgen
{

    Meshing2 :: Meshing2 (const MeshingParameters & mp, const Box<3> & aboundingbox)
    {
        boundingbox = aboundingbox;
        LoadRules (NULL);

        adfront = new AdFront2(boundingbox);
        maxarea = -1;
    }


    Meshing2 :: ~Meshing2 ()
    {
        delete adfront;
        for (int i = 0; i < rules.Size(); i++)
            delete rules[i];
    }

    void Meshing2 :: AddPoint (const Point3d & p, PointIndex globind, 
        MultiPointGeomInfo * mgi,
        bool pointonsurface)
    {
        adfront ->AddPoint (p, globind, mgi, pointonsurface);
    }

    void Meshing2 :: AddBoundaryElement (int i1, int i2,
        const PointGeomInfo & gi1, const PointGeomInfo & gi2)
    {
        if (!gi1.trignum || !gi2.trignum)
        {
            cout << "addboundaryelement: illegal geominfo" << endl;
        }
        adfront -> AddLine (i1-1, i2-1, gi1, gi2);
    }

    void Meshing2 :: StartMesh ()
    {
        foundmap.SetSize (rules.Size());
        canuse.SetSize (rules.Size());
        ruleused.SetSize (rules.Size());

        foundmap = 0;
        canuse = 0;
        ruleused = 0;
    }

    void Meshing2 :: SetMaxArea (double amaxarea)
    {
        maxarea = amaxarea;
    }

    double Meshing2 :: CalcLocalH (const Point3d & /* p */, double gh) const
    {
        return gh;
    }

    void Meshing2 :: DefineTransformation (const Point3d & p1, const Point3d & p2,
        const PointGeomInfo * geominfo1,
        const PointGeomInfo * geominfo2)
    {
        globp1 = p1;
        ex = p2 - p1;
        ex /= ex.Length();
        ey.X() = -ex.Y();
        ey.Y() =  ex.X();
        ey.Z() = 0;
    }

    void Meshing2 :: TransformToPlain (const Point3d & locpoint, 
        const MultiPointGeomInfo & geominf,
        Point2d & plainpoint, double h, int & zone)
    {
        Vec3d p1p (globp1, locpoint);

        //    p1p = locpoint - globp1;
        p1p /= h;
        plainpoint.X() = p1p * ex;
        plainpoint.Y() = p1p * ey;
        zone = 0;
    }

    int Meshing2 :: TransformFromPlain (Point2d & plainpoint,
        Point3d & locpoint, 
        PointGeomInfo & gi, 
        double h)
    {
        Vec3d p1p;
        gi.trignum = 1;

        p1p = plainpoint.X() * ex + plainpoint.Y() * ey;
        p1p *= h;
        locpoint = globp1 + p1p;
        return 0;
    }

    int Meshing2 :: BelongsToActiveChart (const Point3d & p, 
        const PointGeomInfo & gi)
    {
        return 1;
    }

    int Meshing2 :: ComputePointGeomInfo (const Point3d & p, PointGeomInfo & gi)
    {
        gi.trignum = 1;
        return 0;
    }

    int Meshing2 :: ChooseChartPointGeomInfo (const MultiPointGeomInfo & mpgi, 
        PointGeomInfo & pgi)
    {
        pgi = mpgi.GetPGI(1);
        return 0;
    }

    int Meshing2 :: 
        IsLineVertexOnChart (const Point3d & p1, const Point3d & p2,
        int endpoint, const PointGeomInfo & geominfo)
    {
        return 1;
    }

    void Meshing2 ::
        GetChartBoundary (Array<Point2d> & points, 
        Array<Point3d> & points3d, 
        Array<INDEX_2> & lines, double h) const
    {
        points.SetSize (0);
        points3d.SetSize (0);
        lines.SetSize (0);
    }

    double Meshing2 :: Area () const
    {
        return -1;
    }

    MESHING2_RESULT Meshing2 :: GenerateMesh (Mesh & mesh, const MeshingParameters & mp, double gh, int facenr)
    {

        Array<int> pindex, lindex;
        Array<int> delpoints, dellines;

        Array<PointGeomInfo> upgeominfo;  // unique info
        Array<MultiPointGeomInfo> mpgeominfo;  // multiple info

        Array<Element2d> locelements;

        int z1, z2, oldnp(-1);
        bool found;
        int rulenr(-1);
        Point<3> p1, p2;

        const PointGeomInfo * blgeominfo1;
        const PointGeomInfo * blgeominfo2;

        double h, his, hshould;


        Array<Point3d> locpoints;
        Array<int> legalpoints;
        Array<Point2d> plainpoints;
        Array<int> plainzones;
        Array<INDEX_2> loclines;
        int cntelem = 0, nfaces = 0;
        int oldnl = 0;
        int qualclass;

        // test for 3d overlaps
        Box3dTree surfeltree (boundingbox.PMin(),
            boundingbox.PMax());

        Array<int> intersecttrias;
        Array<Point3d> critpoints;

        // test for doubled edges
        StartMesh();

        Array<Point2d> chartboundpoints;
        Array<Point3d> chartboundpoints3d;
        Array<INDEX_2> chartboundlines;

        // illegal points: points with more then 50 elements per node
        int maxlegalpoint(-1), maxlegalline(-1);
        Array<int,PointIndex::BASE> trigsonnode;
        Array<int,PointIndex::BASE> illegalpoint;

        trigsonnode.SetSize (mesh.GetNP());
        illegalpoint.SetSize (mesh.GetNP());

        trigsonnode = 0;
        illegalpoint = 0;

        double totalarea = Area ();
        double meshedarea = 0;

        Array<SurfaceElementIndex> seia;
        mesh.GetSurfaceElementsOfFace (facenr, seia);
        for (int i = 0; i < seia.Size(); i++)
        {
            const Element2d & sel = mesh[seia[i]];

            if (sel.IsDeleted()) continue;

            Box<3> box;
            box.Set ( mesh[sel[0]] );
            box.Add ( mesh[sel[1]] );
            box.Add ( mesh[sel[2]] );
            surfeltree.Insert (box, seia[i]);
        }

        if (totalarea > 0 || maxarea > 0)
            meshedarea = mesh.SurfaceArea();

        adfront ->SetStartFront ();

        double meshedarea_before = meshedarea;

        while (!adfront ->Empty())
        {
            // known for STL meshing
            locpoints.SetSize(0);
            loclines.SetSize(0);
            pindex.SetSize(0);
            lindex.SetSize(0);
            delpoints.SetSize(0);
            dellines.SetSize(0);
            locelements.SetSize(0);

            // unique-pgi, multi-pgi
            upgeominfo.SetSize(0);
            mpgeominfo.SetSize(0);

            nfaces = adfront->GetNFL();

            int baselineindex = adfront -> SelectBaseLine (p1, p2, blgeominfo1, blgeominfo2, qualclass);

            found = 1;

            his = Dist (p1, p2);

            Point3d pmid = Center (p1, p2);
            hshould = CalcLocalH (pmid, mesh.GetH (pmid));
            if (gh < hshould) hshould = gh;

            mesh.RestrictLocalH (pmid, hshould);

            h = hshould;

            double hinner = (3 + qualclass) * max2 (his, hshould);

            adfront ->GetLocals (baselineindex, locpoints, mpgeominfo, loclines, 
                pindex, lindex, 2*hinner);

            if (qualclass > mp.giveuptol2d)
            {
                break;
            }

            PointIndex gpi1 = adfront -> GetGlobalIndex (pindex.Get(loclines[0].I1()));
            PointIndex gpi2 = adfront -> GetGlobalIndex (pindex.Get(loclines[0].I2()));

            // problem recognition !
            if (found && 
                (gpi1 < illegalpoint.Size()+PointIndex::BASE) && 
                (gpi2 < illegalpoint.Size()+PointIndex::BASE) )
            {
                if (illegalpoint[gpi1] || illegalpoint[gpi2])
                    found = 0;
            }

            Point2d p12d, p22d;

            if (found)
            {
                oldnp = locpoints.Size();
                oldnl = loclines.Size();

                DefineTransformation (p1, p2, blgeominfo1, blgeominfo2);

                plainpoints.SetSize (locpoints.Size());
                plainzones.SetSize (locpoints.Size());

                for (int i = 1; i <= locpoints.Size(); i++)
                {
                    TransformToPlain (locpoints.Get(i), 
                        mpgeominfo.Get(i),
                        plainpoints.Elem(i), h, plainzones.Elem(i));

                }

                p12d = plainpoints.Get(1);
                p22d = plainpoints.Get(2);

                for (int i = 2; i <= loclines.Size(); i++)  // don't remove first line
                {
                    z1 = plainzones.Get(loclines.Get(i).I1());
                    z2 = plainzones.Get(loclines.Get(i).I2());

                    // one inner point, one outer
                    if ( (z1 >= 0) != (z2 >= 0))
                    {
                        int innerp = (z1 >= 0) ? 1 : 2;
                        if (IsLineVertexOnChart (locpoints.Get(loclines.Get(i).I1()),
                            locpoints.Get(loclines.Get(i).I2()),
                            innerp,
                            adfront->GetLineGeomInfo (lindex.Get(i), innerp)))
                            // pgeominfo.Get(loclines.Get(i).I(innerp))))
                        {		
                            // use one end of line
                            int pini, pouti;
                            Vec2d v;

                            pini = loclines.Get(i).I(innerp);
                            pouti = loclines.Get(i).I(3-innerp);

                            Point2d pin (plainpoints.Get(pini));
                            Point2d pout (plainpoints.Get(pouti));
                            v = pout - pin;
                            double len = v.Length();
                            if (len <= 1e-6)
                                cout << "WARNING(js): inner-outer: short vector" << endl;
                            else
                                v /= len;

                            Point2d newpout = pin + 1000 * v;
                            newpout = pout;


                            plainpoints.Append (newpout);
                            Point3d pout3d = locpoints.Get(pouti);
                            locpoints.Append (pout3d);

                            plainzones.Append (0);
                            pindex.Append (-1);
                            oldnp++;
                            loclines.Elem(i).I(3-innerp) = oldnp;

                        }
                        else
                        {
                            // remove line
                            loclines.DeleteElement(i);
                            lindex.DeleteElement(i);
                            oldnl--;
                            i--;
                        }			
                    }

                    else if ( (z1 > 0 && z2 > 0 && (z1 != z2)) || ((z1 < 0) && (z2 < 0)) )
                    {
                        loclines.DeleteElement(i);
                        lindex.DeleteElement(i);
                        oldnl--;
                        i--;
                    }
                }

                legalpoints.SetSize(plainpoints.Size());
                for (int i = 1; i <= legalpoints.Size(); i++)
                    legalpoints.Elem(i) = 1;

                double avy = 0;
                for (int i = 1; i <= plainpoints.Size(); i++)
                    avy += plainpoints.Elem(i).Y();
                avy *= 1.0/plainpoints.Size();

                for (int i = 1; i <= plainpoints.Size(); i++)
                {
                    if (plainzones.Elem(i) < 0)
                    {
                        plainpoints.Elem(i) = Point2d (1e4, 1e4);
                        legalpoints.Elem(i) = 0;
                    }
                    if (pindex.Elem(i) == -1)
                    {
                        legalpoints.Elem(i) = 0;
                    }


                    if (plainpoints.Elem(i).Y() < -1e-10*avy) // changed
                    {
                        legalpoints.Elem(i) = 0;
                    }
                }

                GetChartBoundary (chartboundpoints, 
                    chartboundpoints3d,
                    chartboundlines, h);

                oldnp = plainpoints.Size();

                maxlegalpoint = locpoints.Size();
                maxlegalline = loclines.Size();

                if (mp.checkchartboundary)
                {
                    for (int i = 1; i <= chartboundpoints.Size(); i++)
                    {
                        plainpoints.Append (chartboundpoints.Get(i));
                        locpoints.Append (chartboundpoints3d.Get(i));
                        legalpoints.Append (0);
                    }

                    for (int i = 1; i <= chartboundlines.Size(); i++)
                    {
                        INDEX_2 line (chartboundlines.Get(i).I1()+oldnp,
                            chartboundlines.Get(i).I2()+oldnp);
                        loclines.Append (line);
                    }
                }

                oldnl = loclines.Size();
                oldnp = plainpoints.Size();
            }


            if (found)
            {
                rulenr = ApplyRules (plainpoints, legalpoints, maxlegalpoint,
                    loclines, maxlegalline, locelements,
                    dellines, qualclass, mp);

                if (!rulenr)
                {
                    found = 0;
                }
            }

            for (int i = 1; i <= locelements.Size() && found; i++)
            {
                const Element2d & el = locelements.Get(i);

                for (int j = 1; j <= el.GetNP(); j++)
                    if (el.PNum(j) <= oldnp && pindex.Get(el.PNum(j)) == -1)
                    {
                        found = 0;
                        cout << "meshing2, index missing" << endl;
                    }
            }

            if (found)
            {
                locpoints.SetSize (plainpoints.Size());
                upgeominfo.SetSize(locpoints.Size());

                for (int i = oldnp+1; i <= plainpoints.Size(); i++)
                {
                    int err =
                        TransformFromPlain (plainpoints.Elem(i), locpoints.Elem(i), 
                        upgeominfo.Elem(i), h);

                    if (err)
                    {
                        found = 0;
                        break;
                    }
                }
            }

            if (found) 
            {
                double violateminh = 3 + 0.1 * sqr (qualclass);
                double minh = 1e8;
                double newedgemaxh = 0;
                for (int i = oldnl+1; i <= loclines.Size(); i++)
                {
                    double eh = Dist (locpoints.Get(loclines.Get(i).I1()),
                        locpoints.Get(loclines.Get(i).I2()));

                    if (eh > newedgemaxh)
                        newedgemaxh = eh;
                }

                for (int i = 1; i <= locelements.Size(); i++)
                {
                    Point3d pmin = locpoints.Get(locelements.Get(i).PNum(1));
                    Point3d pmax = pmin;
                    for (int j = 2; j <= locelements.Get(i).GetNP(); j++)
                    {
                        const Point3d & hp = 
                            locpoints.Get(locelements.Get(i).PNum(j));
                        pmin.SetToMin (hp);
                        pmax.SetToMax (hp);
                    }
                    double eh = mesh.GetMinH (pmin, pmax);
                    if (eh < minh)
                        minh = eh;
                }

                for (int i = 1; i <= locelements.Size(); i++)
                    for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
                        if (Dist2 (locpoints.Get(locelements.Get(i).PNum(j)), pmid) > hinner*hinner)
                            found = 0;

                static double maxviolate = 0;
                if (newedgemaxh / minh > maxviolate)
                {
                    maxviolate = newedgemaxh / minh;
                }

                if (newedgemaxh > violateminh * minh)
                {
                    found = 0;
                    loclines.SetSize (oldnl);
                    locpoints.SetSize (oldnp);
                }
            }

            // changed for OCC meshing
            if (found)
            {
                for (int i = 1; i <= locelements.Size(); i++)
                    for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
                    {
                        int pi = locelements.Get(i).PNum(j);
                        if (pi <= oldnp)
                        {

                            if (ChooseChartPointGeomInfo (mpgeominfo.Get(pi), upgeominfo.Elem(pi)))
                            {
                                // cannot select, compute new one
                                if (ComputePointGeomInfo (locpoints.Get(pi), upgeominfo.Elem(pi)))
                                {
                                    found = 0;
                                    cout <<("meshing2d, geominfo failed");
                                }
                            }
                        }
                    }
            }

            if (found && mp.checkoverlap)
            {
                // test for overlaps

                Point3d hullmin(1e10, 1e10, 1e10);
                Point3d hullmax(-1e10, -1e10, -1e10);

                for (int i = 1; i <= locelements.Size(); i++)
                    for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
                    {
                        const Point3d & p = locpoints.Get(locelements.Get(i).PNum(j));
                        hullmin.SetToMin (p);
                        hullmax.SetToMax (p);
                    }
                    hullmin += Vec3d (-his, -his, -his);
                    hullmax += Vec3d ( his,  his,  his);

                    surfeltree.GetIntersecting (hullmin, hullmax, intersecttrias);

                    critpoints.SetSize (0);
                    for (int i = oldnp+1; i <= locpoints.Size(); i++)
                        critpoints.Append (locpoints.Get(i));

                    for (int i = 1; i <= locelements.Size(); i++)
                    {
                        const Element2d & tri = locelements.Get(i);
                        if (tri.GetNP() == 3)
                        {
                            const Point3d & tp1 = locpoints.Get(tri.PNum(1));
                            const Point3d & tp2 = locpoints.Get(tri.PNum(2));
                            const Point3d & tp3 = locpoints.Get(tri.PNum(3));

                            Vec3d tv1 (tp1, tp2);
                            Vec3d tv2 (tp1, tp3);

                            double lam1, lam2;
                            for (lam1 = 0.2; lam1 <= 0.8; lam1 += 0.2)
                                for (lam2 = 0.2; lam2 + lam1 <= 0.8; lam2 += 0.2)
                                {
                                    Point3d hp = tp1 + lam1 * tv1 + lam2 * tv2;
                                    critpoints.Append (hp);
                                }
                        }
                        else if (tri.GetNP() == 4)
                        {
                            const Point3d & tp1 = locpoints.Get(tri.PNum(1));
                            const Point3d & tp2 = locpoints.Get(tri.PNum(2));
                            const Point3d & tp3 = locpoints.Get(tri.PNum(3));
                            const Point3d & tp4 = locpoints.Get(tri.PNum(4));

                            double l1, l2;
                            for (l1 = 0.1; l1 <= 0.9; l1 += 0.1)
                                for (l2 = 0.1; l2 <= 0.9; l2 += 0.1)
                                {
                                    Point3d hp;
                                    hp.X() = 
                                        (1-l1)*(1-l2) * tp1.X() +
                                        l1*(1-l2) * tp2.X() +
                                        l1*l2 * tp3.X() +
                                        (1-l1)*l2 * tp4.X();
                                    hp.Y() = 
                                        (1-l1)*(1-l2) * tp1.Y() +
                                        l1*(1-l2) * tp2.Y() +
                                        l1*l2 * tp3.Y() +
                                        (1-l1)*l2 * tp4.Y();
                                    hp.Z() = 
                                        (1-l1)*(1-l2) * tp1.Z() +
                                        l1*(1-l2) * tp2.Z() +
                                        l1*l2 * tp3.Z() +
                                        (1-l1)*l2 * tp4.Z();


                                    critpoints.Append (hp);
                                }
                        }
                    }

                    for (int i = 1; i <= critpoints.Size(); i++)
                    {
                        const Point3d & p = critpoints.Get(i);

                        for (int jj = 0; jj < intersecttrias.Size(); jj++)
                        {
                            SurfaceElementIndex j = intersecttrias[jj];
                            const Element2d & el = mesh[j];

                            int ntrig = (el.GetNP() == 3) ? 1 : 2;

                            int jl;
                            for (jl = 1; jl <= ntrig; jl++)
                            {
                                Point3d tp1, tp2, tp3;

                                if (jl == 1)
                                {
                                    tp1 = mesh.Point(el.PNum(1));
                                    tp2 = mesh.Point(el.PNum(2));
                                    tp3 = mesh.Point(el.PNum(3));
                                }
                                else
                                {
                                    tp1 = mesh.Point(el.PNum(1));
                                    tp2 = mesh.Point(el.PNum(3));
                                    tp3 = mesh.Point(el.PNum(4));
                                }

                                int onchart = 0;
                                for (int k = 1; k <= el.GetNP(); k++)
                                    if (BelongsToActiveChart (mesh.Point(el.PNum(k)),
                                        el.GeomInfoPi(k)))
                                        onchart = 1;
                                if (!onchart)
                                    continue;

                                Vec3d e1(tp1, tp2);
                                Vec3d e2(tp1, tp3);
                                Vec3d n = Cross (e1, e2);
                                n /= n.Length();
                                double lam1, lam2, lam3;
                                lam3 = n * Vec3d (tp1, p);
                                LocalCoordinates (e1, e2, Vec3d (tp1, p), lam1, lam2);

                                if (fabs (lam3) < 0.1 * hshould && 
                                    lam1 > 0 && lam2 > 0 && (lam1 + lam2) < 1)
                                {
                                    for (int k = 1; k <= 5; k++)
                                        adfront -> IncrementClass (lindex.Get(1));

                                    found = 0;

                                }
                            }
                        }
                    }
            }

            if (found)
            {
                // check, whether new front line already exists

                for (int i = oldnl+1; i <= loclines.Size(); i++)
                {
                    int nlgpi1 = loclines.Get(i).I1();
                    int nlgpi2 = loclines.Get(i).I2();
                    if (nlgpi1 <= pindex.Size() && nlgpi2 <= pindex.Size())
                    {
                        nlgpi1 = adfront->GetGlobalIndex (pindex.Get(nlgpi1));
                        nlgpi2 = adfront->GetGlobalIndex (pindex.Get(nlgpi2));

                        int exval = adfront->ExistsLine (nlgpi1, nlgpi2);
                        if (exval)
                        {
                            cout << "ERROR: new line exits, val = " << exval << endl;
                            found = 0;

                        }
                    }
                }
            }

            if (found)
            {
                // everything is ok, perform mesh update
                ruleused.Elem(rulenr)++;

                pindex.SetSize(locpoints.Size());

                for (int i = oldnp+1; i <= locpoints.Size(); i++)
                {
                    PointIndex globind = mesh.AddPoint (locpoints.Get(i));
                    pindex.Elem(i) = adfront -> AddPoint (locpoints.Get(i), globind);
                }

                for (int i = oldnl+1; i <= loclines.Size(); i++)
                {

                    if (pindex.Get(loclines.Get(i).I1()) == -1 || 
                        pindex.Get(loclines.Get(i).I2()) == -1)
                    {
                        cout << "pindex is 0" << endl;
                    }

                    if (!upgeominfo.Get(loclines.Get(i).I1()).trignum || 
                        !upgeominfo.Get(loclines.Get(i).I2()).trignum)
                    {
                        cout << "new el: illegal geominfo" << endl;
                    }

                    adfront -> AddLine (pindex.Get(loclines.Get(i).I1()),
                        pindex.Get(loclines.Get(i).I2()),
                        upgeominfo.Get(loclines.Get(i).I1()),
                        upgeominfo.Get(loclines.Get(i).I2()));
                }
                for (int i = 1; i <= locelements.Size(); i++)
                {
                    Element2d mtri(locelements.Get(i).GetNP());
                    mtri = locelements.Get(i);
                    mtri.SetIndex (facenr);


                    // compute triangle geominfo:
                    for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
                    {
                        mtri.GeomInfoPi(j) = upgeominfo.Get(locelements.Get(i).PNum(j));
                    }

                    for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
                    {
                        mtri.PNum(j) = 
                            locelements.Elem(i).PNum(j) =
                            adfront -> GetGlobalIndex (pindex.Get(locelements.Get(i).PNum(j)));
                    }

                    mesh.AddSurfaceElement (mtri);
                    cntelem++;

                    Box<3> box;
                    box.Set (mesh[mtri[0]]);
                    box.Add (mesh[mtri[1]]);
                    box.Add (mesh[mtri[2]]);
                    surfeltree.Insert (box, mesh.GetNSE()-1);

                    const Point3d & sep1 = mesh.Point (mtri.PNum(1));
                    const Point3d & sep2 = mesh.Point (mtri.PNum(2));
                    const Point3d & sep3 = mesh.Point (mtri.PNum(3));

                    double trigarea = Cross (Vec3d (sep1, sep2), 
                        Vec3d (sep1, sep3)).Length() / 2;

                    if (mtri.GetNP() == 4)
                    {
                        const Point3d & sep4 = mesh.Point (mtri.PNum(4));
                        trigarea += Cross (Vec3d (sep1, sep3), 
                            Vec3d (sep1, sep4)).Length() / 2;
                    }

                    meshedarea += trigarea;

                    if(maxarea > 0 && meshedarea-meshedarea_before > maxarea)
                    {
                        cerr << "meshed area = " << meshedarea-meshedarea_before << endl
                            << "maximal area = " << maxarea << endl
                            << "GIVING UP" << endl;
                        return MESHING2_GIVEUP;
                    }

                    for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
                    {
                        int gpi = locelements.Get(i).PNum(j);

                        int oldts = trigsonnode.Size();
                        if (gpi >= oldts+PointIndex::BASE)
                        {
                            trigsonnode.SetSize (gpi+1-PointIndex::BASE);
                            illegalpoint.SetSize (gpi+1-PointIndex::BASE);
                            for (int k = oldts+PointIndex::BASE; 
                                k <= gpi; k++)
                            {
                                trigsonnode[k] = 0;
                                illegalpoint[k] = 0;
                            }
                        }

                        trigsonnode[gpi]++;

                        if (trigsonnode[gpi] > 20)
                        {
                            illegalpoint[gpi] = 1;
                            cout << "illegal point: " << gpi << endl;
                        }

                        static int mtonnode = 0;
                        if (trigsonnode[gpi] > mtonnode)
                            mtonnode = trigsonnode[gpi];
                    }

                }

                for (int i = 1; i <= dellines.Size(); i++)
                    adfront -> DeleteLine (lindex.Get(dellines.Get(i)));

            }
            else
            {
                adfront -> IncrementClass (lindex.Get(1));

            }

        }

        if (!adfront->Empty())
            return MESHING2_GIVEUP;

        return MESHING2_OK;
    }
}
