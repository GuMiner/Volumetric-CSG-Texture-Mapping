#include "stdafx.h"
#include "mystdlib.h"

#include "myadt.hpp"

#include "linalg.hpp"
#include "csg.hpp"
#include "meshing.hpp"


namespace netgen
{
    Array<SpecialPoint> specpoints;
    static Array<MeshPoint> spoints;

#define TCL_OK 0
#define TCL_ERROR 1

    static void FindPoints (CSGeometry & geom, Mesh & mesh)
    {
        for (int i = 0; i < geom.GetNUserPoints(); i++)
        {
            mesh.AddPoint(geom.GetUserPoint (i));
            mesh.Points().Last().Singularity (geom.GetUserPointRefFactor(i));
            mesh.AddLockedPoint (PointIndex (i+1));
        }

        SpecialPointCalculation spc;

        spc.SetIdEps(geom.GetIdEps());

        if (spoints.Size() == 0)
            spc.CalcSpecialPoints (geom, spoints);

        spc.AnalyzeSpecialPoints (geom, spoints, specpoints);
        if (mparam.enableOutput)
            cout << " " << specpoints.Size() << " special points." << endl;
    }

    static void FindEdges (CSGeometry & geom, Mesh & mesh, const bool setmeshsize = false)
    {
        EdgeCalculation ec (geom, specpoints);
        ec.SetIdEps(geom.GetIdEps());
        ec.Calc (mparam.maxh, mesh);

        for (int i = 0; i < geom.singedges.Size(); i++)
        {
            geom.singedges[i]->FindPointsOnEdge (mesh);
            if(setmeshsize)
                geom.singedges[i]->SetMeshSize(mesh,10.0*geom.BoundingBox().Diam());
        }
        for (int i = 0; i < geom.singpoints.Size(); i++)
            geom.singpoints[i]->FindPoints (mesh);

        for (int i = 1; i <= mesh.GetNSeg(); i++)
        {
            int ok = 0;
            for (int k = 1; k <= mesh.GetNFD(); k++)
                if (mesh.GetFaceDescriptor(k).SegmentFits (mesh.LineSegment(i)))
                {
                    ok = k;
                }

            if (!ok)
            {
                ok = mesh.AddFaceDescriptor (FaceDescriptor (mesh.LineSegment(i)));
            }

            mesh.LineSegment(i).si = ok;
        }

        if (geom.identifications.Size())
        {
            for (int i = 0; i < geom.identifications.Size(); i++)
            {
                geom.identifications[i]->IdentifyPoints (mesh);

            }
            for (int i = 0; i < geom.identifications.Size(); i++)
                geom.identifications[i]->IdentifyFaces (mesh);
        }

        // find intersecting segments
        Point3d pmin, pmax;
        mesh.GetBox (pmin, pmax);
        Box3dTree segtree (pmin, pmax);

        for (SegmentIndex si = 0; si < mesh.GetNSeg(); si++)
        {
            if (mesh[si].seginfo)
            {
                Box<3> hbox;
                hbox.Set (mesh[mesh[si][0]]);
                hbox.Add (mesh[mesh[si][1]]);
                segtree.Insert (hbox.PMin(), hbox.PMax(), si);
            }
        }

        Array<int> loc;
        if (!ec.point_on_edge_problem)
            for (SegmentIndex si = 0; si < mesh.GetNSeg(); si++)
            {
                if (!mesh[si].seginfo) continue;

                Box<3> hbox;
                hbox.Set (mesh[mesh[si][0]]);
                hbox.Add (mesh[mesh[si][1]]);
                hbox.Increase (1e-6);
                segtree.GetIntersecting (hbox.PMin(), hbox.PMax(), loc);

                for (int j = 0; j < loc.Size(); j++)
                {
                    SegmentIndex sj = loc[j];
                    if (sj >= si) continue;
                    if (!mesh[si].seginfo || !mesh[sj].seginfo) continue;
                    if (mesh[mesh[si][0]].GetLayer() != mesh[mesh[sj][1]].GetLayer()) continue;

                    Point<3> pi1 = mesh[mesh[si][0]];
                    Point<3> pi2 = mesh[mesh[si][1]];
                    Point<3> pj1 = mesh[mesh[sj][0]];
                    Point<3> pj2 = mesh[mesh[sj][1]];
                    Vec<3> vi = pi2 - pi1;
                    Vec<3> vj = pj2 - pj1;

                    if (sqr (vi * vj) > (1.-1e-6) * Abs2 (vi) * Abs2 (vj)) continue;

                    // pi1 + vi t = pj1 + vj s
                    Mat<3,2> mat;
                    Vec<3> rhs;
                    Vec<2> sol;

                    for (int jj = 0; jj < 3; jj++)
                    { 
                        mat(jj,0) = vi(jj); 
                        mat(jj,1) = -vj(jj); 
                        rhs(jj) = pj1(jj)-pi1(jj); 
                    }

                    mat.Solve (rhs, sol);

                    if (sol(0) > 1e-6 && sol(0) < 1-1e-6 &&
                        sol(1) > 1e-6 && sol(1) < 1-1e-6 &&
                        Abs (rhs - mat*sol) < 1e-6)
                    {
                        Point<3> ip = pi1 + sol(0) * vi;

                        Point<3> pip = ip;
                        ProjectToEdge (geom.GetSurface (mesh[si].surfnr1),
                            geom.GetSurface (mesh[si].surfnr2), pip);

                        if (Dist (ip, pip) > 1e-6*geom.MaxSize()) continue;
                        pip = ip;
                        ProjectToEdge (geom.GetSurface (mesh[sj].surfnr1),
                            geom.GetSurface (mesh[sj].surfnr2), pip);

                        if (Dist (ip, pip) > 1e-6*geom.MaxSize()) continue;

                        geom.AddUserPoint (ip);
                        spoints.Append (MeshPoint (ip, mesh[mesh[si][0]].GetLayer()));
                        mesh.AddPoint (ip);
                    }
                }
            }  
    }

    static void MeshSurface (CSGeometry & geom, Mesh & mesh)
    {
        Array<Segment> segments;
        int noldp = mesh.GetNP();

        // find master faces from identified
        Array<int> masterface(mesh.GetNFD());
        for (int i = 1; i <= mesh.GetNFD(); i++)
            masterface.Elem(i) = i;

        Array<INDEX_2> fpairs;
        bool changed;
        do
        {
            changed = 0;
            for (int i = 0; i < geom.identifications.Size(); i++)
            {
                geom.identifications[i]->GetIdentifiedFaces (fpairs);

                for (int j = 0; j < fpairs.Size(); j++)
                {
                    if (masterface.Get(fpairs[j].I1()) <
                        masterface.Get(fpairs[j].I2()))
                    {
                        changed = 1;
                        masterface.Elem(fpairs[j].I2()) =
                            masterface.Elem(fpairs[j].I1());
                    }
                    if (masterface.Get(fpairs[j].I2()) <
                        masterface.Get(fpairs[j].I1()))
                    {
                        changed = 1;
                        masterface.Elem(fpairs[j].I1()) =
                            masterface.Elem(fpairs[j].I2());
                    }
                }
            }
        }
        while (changed);


        int bccnt=0;
        for (int k = 0; k < geom.GetNSurf(); k++)
            bccnt = max2 (bccnt, geom.GetSurface(k)->GetBCProperty());

        for (int k = 1; k <= mesh.GetNFD(); k++)
        {
            bool increased = false;

            FaceDescriptor & fd = mesh.GetFaceDescriptor(k);
            const Surface * surf = geom.GetSurface(fd.SurfNr());

            if (fd.TLOSurface() && 
                geom.GetTopLevelObject(fd.TLOSurface()-1) -> GetBCProp() > 0)
                fd.SetBCProperty (geom.GetTopLevelObject(fd.TLOSurface()-1) -> GetBCProp());
            else if (surf -> GetBCProperty() != -1)
                fd.SetBCProperty (surf->GetBCProperty());
            else
            {
                bccnt++;
                fd.SetBCProperty (bccnt);
                increased = true;
            }      

            for (int l = 0; l < geom.bcmodifications.Size(); l++)
            {
                if (geom.GetSurfaceClassRepresentant (fd.SurfNr()) == 
                    geom.GetSurfaceClassRepresentant (geom.bcmodifications[l].si) &&
                    (fd.DomainIn() == geom.bcmodifications[l].tlonr+1 ||
                    fd.DomainOut() == geom.bcmodifications[l].tlonr+1))
                {
                    if(geom.bcmodifications[l].bcname == NULL)
                        fd.SetBCProperty (geom.bcmodifications[l].bcnr);
                    else
                    {
                        if(!increased)
                        {
                            bccnt++;
                            fd.SetBCProperty (bccnt);
                            increased = true;
                        }
                    }
                }
            }
        }

        mesh.SetNBCNames( bccnt );

        for (int k = 1; k <= mesh.GetNFD(); k++)
        {
            FaceDescriptor & fd = mesh.GetFaceDescriptor(k);
            const Surface * surf = geom.GetSurface(fd.SurfNr());
            if (fd.TLOSurface() )
            {
                int bcp = fd.BCProperty();
                string nextbcname = geom.GetTopLevelObject(fd.TLOSurface()-1) -> GetBCName();
                if ( nextbcname != "default" )
                    mesh.SetBCName ( bcp - 1 , nextbcname );
            }
            else // if (surf -> GetBCProperty() != -1)
            {
                int bcp = fd.BCProperty();
                string nextbcname = surf->GetBCName();
                if ( nextbcname != "default" )
                    mesh.SetBCName ( bcp - 1, nextbcname );
            }
        }

        for (int k = 1; k <= mesh.GetNFD(); k++)
        {
            FaceDescriptor & fd = mesh.GetFaceDescriptor(k);
            fd.SetBCName ( mesh.GetBCNamePtr ( fd.BCProperty() - 1 ) );
        }


        //!!

        for (int k = 1; k <= mesh.GetNFD(); k++)
        {
            FaceDescriptor & fd = mesh.GetFaceDescriptor(k);
            //const Surface * surf = geom.GetSurface(fd.SurfNr());

            for (int l = 0; l < geom.bcmodifications.Size(); l++)
            {
                if (geom.GetSurfaceClassRepresentant (fd.SurfNr()) == 
                    geom.GetSurfaceClassRepresentant (geom.bcmodifications[l].si) &&
                    (fd.DomainIn() == geom.bcmodifications[l].tlonr+1 ||
                    fd.DomainOut() == geom.bcmodifications[l].tlonr+1) &&
                    geom.bcmodifications[l].bcname != NULL
                    )
                {
                    int bcp = fd.BCProperty();
                    mesh.SetBCName ( bcp - 1, *(geom.bcmodifications[l].bcname) );
                    fd.SetBCName ( mesh.GetBCNamePtr ( bcp - 1) );
                }
            }
        }

        for(int k = 0; k<geom.bcmodifications.Size(); k++)
        {
            delete geom.bcmodifications[k].bcname;
            geom.bcmodifications[k].bcname = NULL;
        }

        //!!


        for (int j = 0; j < geom.singfaces.Size(); j++)
        {
            Array<int> surfs;
            geom.GetIndependentSurfaceIndices (geom.singfaces[j]->GetSolid(),
                geom.BoundingBox(), surfs);
            for (int k = 1; k <= mesh.GetNFD(); k++)
            {
                FaceDescriptor & fd = mesh.GetFaceDescriptor(k);
                for (int l = 0; l < surfs.Size(); l++)
                    if (surfs[l] == fd.SurfNr())
                    {
                        if (geom.singfaces[j]->GetDomainNr() == fd.DomainIn())
                            fd.SetDomainInSingular (1);
                        if (geom.singfaces[j]->GetDomainNr() == fd.DomainOut())
                            fd.SetDomainOutSingular (1);
                    }
            }
        }


        // assemble edge hash-table
        mesh.CalcSurfacesOfNode();

        for (int k = 1; k <= mesh.GetNFD(); k++)
        {
            if (masterface.Get(k) != k)
                continue;

            FaceDescriptor & fd = mesh.GetFaceDescriptor(k);
            if (mparam.enableOutput)
                cout << "Surface " << k << endl;

            int oldnf = mesh.GetNSE();

            const Surface * surf =
                geom.GetSurface((mesh.GetFaceDescriptor(k).SurfNr()));

            Meshing2Surfaces meshing(*surf, mparam, geom.BoundingBox());

            double eps = 1e-8 * geom.MaxSize();
            for (PointIndex pi = PointIndex::BASE; pi < noldp+PointIndex::BASE; pi++)
            {
                meshing.AddPoint (mesh[pi], pi, NULL,
                    (surf->PointOnSurface(mesh[pi], eps) != 0));
            }

            segments.SetSize (0);

            for (SegmentIndex si = 0; si < mesh.GetNSeg(); si++)
                if (mesh[si].si == k)
                {
                    segments.Append (mesh[si]);
                }

            for (int i = 1; i <= geom.identifications.Size(); i++)
            {
                geom.identifications.Get(i)->
                    BuildSurfaceElements(segments, mesh, surf);
            }

            for (int si = 0; si < segments.Size(); si++)
            {
                PointGeomInfo gi;
                gi.trignum = k;
                meshing.AddBoundaryElement (segments[si][0] + 1 - PointIndex::BASE, 
                    segments[si][1] + 1 - PointIndex::BASE, 
                    gi, gi);
            }

            double maxh = mparam.maxh;
            if (fd.DomainIn() != 0)
            {
                const Solid * s1 = 
                    geom.GetTopLevelObject(fd.DomainIn()-1) -> GetSolid();
                if (s1->GetMaxH() < maxh)
                    maxh = s1->GetMaxH();
                maxh = min2(maxh, geom.GetTopLevelObject(fd.DomainIn()-1)->GetMaxH());
            }
            if (fd.DomainOut() != 0)
            {
                const Solid * s1 = 
                    geom.GetTopLevelObject(fd.DomainOut()-1) -> GetSolid();
                if (s1->GetMaxH() < maxh)
                    maxh = s1->GetMaxH();
                maxh = min2(maxh, geom.GetTopLevelObject(fd.DomainOut()-1)->GetMaxH());
            }
            if (fd.TLOSurface() != 0)
            {
                double hi = geom.GetTopLevelObject(fd.TLOSurface()-1) -> GetMaxH();
                if (hi < maxh) maxh = hi;
            }

            mparam.checkoverlap = 0;

            MESHING2_RESULT res =
                meshing.GenerateMesh (mesh, mparam, maxh, k);

            if (res != MESHING2_OK)
            {
                throw NgException ("Problem in Surface mesh generation");
            }


            for (SurfaceElementIndex sei = oldnf; sei < mesh.GetNSE(); sei++)
                mesh[sei].SetIndex (k);

            if (mparam.enableOutput)
                cout << "  " << (mesh.GetNSE() - oldnf) << " elements, " << mesh.GetNP() << " points" << std::endl;

        }

        mesh.Compress();

        do
        {
            changed = 0;
            for (int k = 1; k <= mesh.GetNFD(); k++)
            {

                if (masterface.Get(k) == k)
                    continue;

                FaceDescriptor & fd = mesh.GetFaceDescriptor(k);
                if (mparam.enableOutput)
                    cout << "Surface " << k << endl;

                int oldnf = mesh.GetNSE();

                const Surface * surf =
                    geom.GetSurface((mesh.GetFaceDescriptor(k).SurfNr()));

                segments.SetSize (0);
                for (int i = 1; i <= mesh.GetNSeg(); i++)
                {
                    Segment * seg = &mesh.LineSegment(i);
                    if (seg->si == k)
                        segments.Append (*seg);
                }

                for (int i = 1; i <= geom.identifications.Size(); i++)
                {
                    geom.identifications.Elem(i)->GetIdentifiedFaces (fpairs);
                    int found = 0;
                    for (int j = 1; j <= fpairs.Size(); j++)
                        if (fpairs.Get(j).I1() == k || fpairs.Get(j).I2() == k)
                            found = 1;

                    if (!found)
                        continue;

                    geom.identifications.Get(i)->
                        BuildSurfaceElements(segments, mesh, surf);
                    if (!segments.Size())
                        break;
                }

                for (SurfaceElementIndex  sei = oldnf; sei < mesh.GetNSE(); sei++)
                    mesh[sei].SetIndex (k);


                if (!segments.Size())
                {
                    masterface.Elem(k) = k;
                    changed = 1; 
                }
                if (mparam.enableOutput)
                    cout << "  " << (mesh.GetNSE() - oldnf) << " elements, " << mesh.GetNP() << " points" << std::endl;
            }

        }
        while (changed);


        mesh.SplitSeparatedFaces();
        mesh.CalcSurfacesOfNode();
    }

    static double TriangleQualityInst (const Point3d & p1, const Point3d & p2,
        const Point3d & p3)
    {
        // quality 0 (worst) .. 1 (optimal)

        Vec3d v1, v2, v3;
        double s1, s2, s3;
        double an1, an2, an3;

        v1 = p2 - p1;
        v2 = p3 - p1;
        v3 = p3 - p2;

        an1 = Angle (v1, v2);
        v1 *= -1;
        an2 = Angle (v1, v3);
        an3 = Angle (v2, v3);

        s1 = sin (an1/2);
        s2 = sin (an2/2);
        s3 = sin (an3/2);

        return 8 * s1 * s2 * s3;
    }

    // Displays mesh quality information
    void MeshQuality2d (const Mesh& mesh)
    {
        if (!mparam.enableOutput)
            return;

        // Number of quality levels and stars (max)
        int ncl = 10, stars = 30;
        
        Array<INDEX> incl(ncl);
        incl = 0;
        for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
        {
            double qual = TriangleQualityInst (mesh[mesh[sei][0]],
                mesh[mesh[sei][1]],
                mesh[mesh[sei][2]]);

            int cl = int ( (ncl-1e-3) * qual ) + 1;
            incl.Elem(cl)++;
        }

        cout << endl;

        cout << "Points:    " << mesh.GetNP() << endl;
        cout << "Triangles: " << mesh.GetNSE() << endl;

        cout << endl;
        cout << "Quality Distribution (%):" << endl;
        
        for (INDEX i = ncl; i > 0; i--)
        {
            cout << " " << i*100/ncl << "-" << (i-1)*100/ncl << ": ";

            for (int j = 0; j < (int)(stars*(double)incl.Get(i)/(double)mesh.GetNSE()); j++)
            {
                cout << "*";
            }
            cout << endl;
        }
        cout << endl;
    }

    /**
    * Generates a surface mesh. If mesh is null, a new mesh is created -- else the current mesh is deleted.
    * Runs an Analyze-Edge Mesh-Surface Mesh procedure.
    **/
    int CSGGenerateMesh (CSGeometry & geom, Mesh *& mesh, MeshingParameters & mparam)
    {
        if (mesh)
        {
            mesh->DeleteMesh();
        }
        else
        {
            mesh = new Mesh();
        }

        mesh->SetGlobalH (mparam.maxh);
        mesh->SetMinimalH (mparam.minh);

        Array<double> maxhdom(geom.GetNTopLevelObjects());
        for (int i = 0; i < maxhdom.Size(); i++)
        {
            maxhdom[i] = geom.GetTopLevelObject(i)->GetMaxH();
        }

        mesh->SetMaxHDomain (maxhdom);
        if (mparam.uselocalh)
        {
            double maxsize = geom.MaxSize(); 
            mesh->SetLocalH (Point<3>(-maxsize, -maxsize, -maxsize),
                Point<3>(maxsize, maxsize, maxsize),
                mparam.grading);
        }

        spoints.SetSize(0);

        FindPoints (geom, *mesh);
        if (mparam.enableOutput)
            cout << endl << "Find points done." << endl;


        FindEdges (geom, *mesh, true);
        if (mparam.enableOutput)
            cout << endl << "Find edges done." << endl;

        if (mparam.uselocalh)
        {
            mesh->CalcLocalH(mparam.grading);
            mesh->DeleteMesh();

            FindPoints (geom, *mesh);
            FindEdges (geom, *mesh, true);

            mesh->DeleteMesh();

            FindPoints (geom, *mesh);
            FindEdges (geom, *mesh);
        }



        MeshSurface (geom, *mesh);
        if (mparam.enableOutput)
            cout << endl << "Mesh surface done." << endl;


        if (mparam.uselocalh && 0)
        {
            mesh->CalcLocalH(mparam.grading);
            mesh->DeleteMesh();

            FindPoints (geom, *mesh);
            FindEdges (geom, *mesh);

            MeshSurface (geom, *mesh);
        }


        MeshQuality2d (*mesh);
        mesh->CalcSurfacesOfNode();

        return TCL_OK;
    }
}
