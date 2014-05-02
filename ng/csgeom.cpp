#include "stdafx.h"
#include "mystdlib.h"
#include "myadt.hpp"

#include "linalg.hpp"
#include "csg.hpp"

namespace netgen
{
    // Meshing parameters for global access
    MeshingParameters mparam;

    TopLevelObject ::  
        TopLevelObject (Solid * asolid,
        Surface * asurface)
    {
        solid = asolid;
        surface = asurface;

        SetLayer (1);

        if (!surface)
            maxh = solid->GetMaxH();
        else
            maxh = surface->GetMaxH();

        SetBCProp (-1);

        bcname = "default";
    }

    Box<3> CSGeometry::default_boundingbox (Point<3> (-1000, -1000, -1000),
        Point<3> ( 1000,  1000,  1000));


    CSGeometry :: CSGeometry ()
        : boundingbox (default_boundingbox),
        identicsurfaces (100), ideps(1e-9), filename("")
    {
        ;
    }

    CSGeometry :: CSGeometry (const string & afilename)
        : boundingbox (default_boundingbox),
        identicsurfaces (100), ideps(1e-9), filename(afilename)
    {
    }

    CSGeometry :: ~CSGeometry ()
    {
        Clean();
    }

    void CSGeometry :: Clean ()
    {
        Array< Solid* > to_delete;

        for (int i = 0; i < solids.Size(); i++)
            if(!to_delete.Contains(solids[i]->S1()))
                to_delete.Append(solids[i]->S1());
        for (int i = 0; i < solids.Size(); i++)
            if(!to_delete.Contains(solids[i]))
                to_delete.Append(solids[i]);

        for(int i = 0; i < to_delete.Size(); i++)
            delete to_delete[i];    

        solids.DeleteAll ();

        for (int i = 0; i < splinecurves2d.Size(); i++)
            delete splinecurves2d[i];
        splinecurves2d.DeleteAll();
        
        for(int i = 0; i<delete_them.Size(); i++)
            delete delete_them[i];
        delete_them.DeleteAll();
        surfaces.DeleteAll();

        for (int i = 0; i < toplevelobjects.Size(); i++)
            delete toplevelobjects[i];
        toplevelobjects.DeleteAll ();

        for(int i = 0; i < identifications.Size(); i++)
            delete identifications[i];
        identifications.DeleteAll();

        for (int i = 0; i < singfaces.Size(); i++)
            delete singfaces[i];
        singfaces.DeleteAll();
        for (int i = 0; i < singedges.Size(); i++)
            delete singedges[i];
        singedges.DeleteAll();
        for (int i = 0; i < singpoints.Size(); i++)
            delete singpoints[i];
        singpoints.DeleteAll();

    }

    extern int CSGGenerateMesh (CSGeometry & geom, Mesh *& mesh, MeshingParameters & mparam);

    int CSGeometry :: GenerateMesh (Mesh*& mesh, MeshingParameters & mparam_loc)
    {
        mparam = mparam_loc;
        return CSGGenerateMesh (*this, mesh, mparam_loc);
    }

    void CSGeometry :: Load (istream & ist)
    {
        char key[100], name[100], classname[100], sname[100];
        int ncoeff, i, j;
        Array<double> coeff;

        while (ist.good())
        {
            ist >> key;
            if (strcmp (key, "boundingbox") == 0)
            {
                Point<3> pmin, pmax;
                ist >> pmin(0) >> pmin(1) >> pmin(2);
                ist >> pmax(0) >> pmax(1) >> pmax(2);
                SetBoundingBox (Box<3> (pmin, pmax));
            }
            if (strcmp (key, "primitive") == 0)
            {
                ist >> name >> classname >> ncoeff;
                coeff.SetSize (ncoeff);
                for (i = 0; i < ncoeff; i++)
                    ist >> coeff[i];

                Primitive * nprim = Primitive::CreatePrimitive (classname);
                nprim -> SetPrimitiveData (coeff);
                Solid * nsol = new Solid (nprim);

                for (j = 0; j < nprim->GetNSurfaces(); j++)
                {
                    sprintf (sname, "%s,%d", name, j);
                    AddSurface (sname, &nprim->GetSurface(j));
                    nprim -> SetSurfaceId (j, GetNSurf());
                }
                SetSolid (name, nsol);
            }
            else if (strcmp (key, "solid") == 0)
            {
                ist >> name;
                Solid * nsol = Solid::CreateSolid (ist, solids);

                cout << " I have found solid " << name << " = ";
                nsol -> GetSolidData (cout);
                cout << endl;

                SetSolid (name, nsol);
            }
            else if (strcmp (key, "toplevel") == 0)
            {
                char type[20], solname[50], surfname[50];
                const Solid * sol = NULL;
                const Surface * surf = NULL;
                int nr;

                ist >> type;
                if (strcmp (type, "solid") == 0)
                {
                    ist >> solname;
                    sol = GetSolid (solname);
                }
                if (strcmp (type, "surface") == 0)
                {
                    ist >> solname >> surfname;
                    sol = GetSolid (solname);
                    surf = GetSurface (surfname);
                }
                nr = SetTopLevelObject ((Solid*)sol, (Surface*)surf);
                //GetTopLevelObject(nr)->SetData(ist);
            }
            else if (strcmp (key, "identify") == 0)
            {
                char type[10], surfname1[50], surfname2[50];
                const Surface * surf1;
                const Surface * surf2;


                ist >> type >> surfname1 >> surfname2;
                surf1 = GetSurface(surfname1);
                surf2 = GetSurface(surfname2);

                AddIdentification (new PeriodicIdentification 
                    (GetNIdentifications(),
                    *this, surf1, surf2));
            }
            else if (strcmp (key, "end") == 0)
                break;
        }
    }

    void CSGeometry :: LoadSurfaces (istream & in)
    {
        Array<double> coeffs;
        string classname;
        int nsurfaces,size;

        in >> classname;

        if(classname == "csgsurfaces")
            in >> nsurfaces;
        else
            nsurfaces = atoi(classname.c_str());

        Point<3> dummypoint(0,0,0);
        Vec<3> dummyvec(0,0,0);
        double dummydouble(0.1);

        for(int i=0; i<nsurfaces; i++)
        {
            in >> classname;
            in >> size;

            coeffs.SetSize(size);

            for(int j=0; j<size; j++)
                in >> coeffs[j];

            if(classname == "plane")
            {
                Plane * plane = new Plane(dummypoint,dummyvec);
                plane->SetPrimitiveData(coeffs);

                AddSurface(plane);
                delete_them.Append(plane);
            }

            else if(classname == "sphere")
            {
                Sphere * sphere = new Sphere(dummypoint,dummydouble);
                sphere->SetPrimitiveData(coeffs);

                AddSurface(sphere);
                delete_them.Append(sphere);
            }

            else if(classname == "cylinder")
            {
                Cylinder * cylinder = new Cylinder(coeffs);

                AddSurface(cylinder);
                delete_them.Append(cylinder);
            }

            else if(classname == "ellipticcylinder")
            {
                EllipticCylinder * cylinder = new EllipticCylinder(coeffs);
                AddSurface(cylinder);
                delete_them.Append(cylinder);
            }


            else if(classname == "torus")
            {
                Torus * torus = new Torus(dummypoint,dummyvec,dummydouble, dummydouble);
                torus->SetPrimitiveData(coeffs);
                AddSurface(torus);
                delete_them.Append(torus);
            }


            else if(classname == "cone")
            {
                Cone * cone = new Cone(dummypoint,dummypoint,dummydouble,dummydouble);
                cone->SetPrimitiveData(coeffs);

                AddSurface(cone);
                delete_them.Append(cone);
            }

            else if(classname == "extrusionface")
            {
                ExtrusionFace * ef =
                    new ExtrusionFace(coeffs);

                AddSurface(ef);
                delete_them.Append(ef);
            }

            else if(classname == "revolutionface")
            {
                RevolutionFace * rf =
                    new RevolutionFace(coeffs);

                AddSurface(rf);
                delete_them.Append(rf);
            }

        }
    }

    void CSGeometry :: AddSurface (Surface * surf)
    {
        static int cntsurfs = 0;
        cntsurfs++;
        char name[15];
        sprintf (name, "nnsurf %d", cntsurfs);
        AddSurface (name, surf);
    }

    void CSGeometry :: AddSurface (char * name, Surface * surf)
    { 
        if (mparam.enableOutput)
            cout << "Adding surface " << name << endl;
        surfaces.Set (name, surf); 
        surf->SetName (name);
    }

    void CSGeometry :: AddSurfaces (Primitive * prim)
    {
        for (int i = 0; i < prim->GetNSurfaces(); i++)
        {
            AddSurface (&prim->GetSurface(i));
            prim->SetSurfaceId (i, GetNSurf()-1);
            surf2prim.Append (prim);
        }
    }

    const Surface * CSGeometry :: GetSurface (const char * name) const
    {
        if (surfaces.Used(name))
            return surfaces.Get(name);
        else
            return NULL;
    }

    void CSGeometry :: SetSolid (const char * name, Solid * sol)
    {
        Solid * oldsol = NULL;

        if (solids.Used (name))
            oldsol = solids.Get(name);

        solids.Set (name, sol);
        sol->SetName (name);

        if (oldsol)
        {
            if (oldsol->op != Solid::ROOT ||
                sol->op != Solid::ROOT)
            {
                cerr << "Setsolid: old or new no root" << endl;
            }
            oldsol -> s1 = sol -> s1;
        }
    }

    const Solid * CSGeometry :: GetSolid (const char * name) const
    {
        if (solids.Used(name))
            return solids.Get(name);
        else
            return NULL;
    }

    const Solid * CSGeometry :: GetSolid (const string & name) const
    {
        if (solids.Used(name.c_str()))
            return solids.Get(name.c_str());
        else
            return NULL;
    }

    void CSGeometry :: SetSplineCurve (const char * name, SplineGeometry<2> * spl)
    {
        splinecurves2d.Set(name,spl);
    }
    void CSGeometry :: SetSplineCurve (const char * name, SplineGeometry<3> * spl)
    {
        splinecurves3d.Set(name,spl);
    }


    const SplineGeometry<2> * CSGeometry :: GetSplineCurve2d (const string & name) const
    {
        if (splinecurves2d.Used(name.c_str()))
            return splinecurves2d.Get(name.c_str());
        else
            return NULL;
    }
    const SplineGeometry<3> * CSGeometry :: GetSplineCurve3d (const string & name) const
    {
        if (splinecurves3d.Used(name.c_str()))
            return splinecurves3d.Get(name.c_str());
        else
            return NULL;
    }

    int CSGeometry :: SetTopLevelObject (Solid * sol, Surface * surf)
    {
        return toplevelobjects.Append (new TopLevelObject (sol, surf)) - 1;
    }

    TopLevelObject * CSGeometry :: 
        GetTopLevelObject (const Solid * sol, const Surface * surf)
    {
        for (int i = 0; i < toplevelobjects.Size(); i++)
        {
            if (toplevelobjects[i]->GetSolid() == sol &&
                toplevelobjects[i]->GetSurface() == surf)
                return (toplevelobjects[i]);
        }
        return NULL;
    }

    void CSGeometry :: RemoveTopLevelObject (Solid * sol, Surface * surf)
    {
        for (int i = 0; i < toplevelobjects.Size(); i++)
        {
            if (toplevelobjects[i]->GetSolid() == sol &&
                toplevelobjects[i]->GetSurface() == surf)
            {
                delete toplevelobjects[i];
                toplevelobjects.DeleteElement (i+1);
                break;
            }
        }
    }

    void CSGeometry :: AddIdentification (Identification * ident)
    {
        identifications.Append (ident);
    }

    void CSGeometry :: SetFlags (const char * solidname, const Flags & flags)
    {
        Solid * solid = solids.Elem(solidname);
        Array<int> surfind;

        int i;
        double maxh = flags.GetNumFlag ("maxh", -1);
        if (maxh > 0 && solid)
        {
            solid->GetSurfaceIndices (surfind);

            for (i = 0; i < surfind.Size(); i++)
            {
                if (surfaces[surfind[i]]->GetMaxH() > maxh)
                    surfaces[surfind[i]] -> SetMaxH (maxh);
            }

            solid->SetMaxH (maxh);
        }

        if ( flags.StringFlagDefined ("bcname") )
        {
            solid->GetSurfaceIndices (surfind);
            string bcn = flags.GetStringFlag("bcname", "default");
            for (i = 0; i < surfind.Size(); i++)
            {
                if(surfaces[surfind[i]]->GetBCName() == "default")
                    surfaces[surfind[i]]->SetBCName(bcn);
            }
        }

        if (flags.StringListFlagDefined ("bcname"))
        {
            const Array<char*> & bcname = flags.GetStringListFlag("bcname");

            Polyhedra * polyh;
            if(solid->S1())
                polyh = dynamic_cast<Polyhedra *>(solid->S1()->GetPrimitive());
            else
                polyh = dynamic_cast<Polyhedra *>(solid->GetPrimitive());

            if(polyh)
            {
                Array < Array<int> * > polysurfs;
                polyh->GetPolySurfs(polysurfs);
                if(bcname.Size() != polysurfs.Size())
                    cerr << "WARNING: solid \"" << solidname << "\" has " << polysurfs.Size()
                    << " surfaces and should get " << bcname.Size() << " bc-names!" << endl;

                for ( i = 0; i < min2(polysurfs.Size(),bcname.Size()); i++)
                {
                    for (int j = 0; j < polysurfs[i]->Size(); j++)
                    {
                        if(surfaces[(*polysurfs[i])[j]]->GetBCName() == "default")
                            surfaces[(*polysurfs[i])[j]]->SetBCName(bcname[i]);
                    }
                    delete polysurfs[i];
                }
            }
            else
            {
                solid->GetSurfaceIndices (surfind);
                if(bcname.Size() != surfind.Size())
                    cerr << "WARNING: solid \"" << solidname << "\" has " << surfind.Size()
                    << " surfaces and should get " << bcname.Size() << " bc-names!" << endl;

                for (i = 0; i < min2(surfind.Size(),bcname.Size()); i++)
                {
                    if(surfaces[surfind[i]]->GetBCName() == "default")
                        surfaces[surfind[i]]->SetBCName(bcname[i]);
                }
            }
        }

        if (flags.NumFlagDefined ("bc"))
        {
            solid->GetSurfaceIndices (surfind);
            int bc = int (flags.GetNumFlag("bc", -1));
            for (i = 0; i < surfind.Size(); i++)
            {
                if (surfaces[surfind[i]]->GetBCProperty() == -1)
                    surfaces[surfind[i]]->SetBCProperty(bc);
            }
        }

        if (flags.NumListFlagDefined ("bc"))
        {
            const Array<double> & bcnum = flags.GetNumListFlag("bc");

            Polyhedra * polyh;
            if(solid->S1())
                polyh = dynamic_cast<Polyhedra *>(solid->S1()->GetPrimitive());
            else
                polyh = dynamic_cast<Polyhedra *>(solid->GetPrimitive());

            if(polyh)
            {
                Array < Array<int> * > polysurfs;
                polyh->GetPolySurfs(polysurfs);
                if(bcnum.Size() != polysurfs.Size())
                    cerr << "WARNING: solid \"" << solidname << "\" has " << polysurfs.Size()
                    << " surfaces and should get " << bcnum.Size() << " bc-numbers!" << endl;

                for ( i = 0; i < min2(polysurfs.Size(),bcnum.Size()); i++)
                {
                    for (int j = 0; j < polysurfs[i]->Size(); j++)
                    {
                        if ( surfaces[(*polysurfs[i])[j]]->GetBCProperty() == -1 )
                            surfaces[(*polysurfs[i])[j]]->SetBCProperty(int(bcnum[i]));
                    }
                    delete polysurfs[i];
                }
            }
            else
            {
                solid->GetSurfaceIndices (surfind);
                if(bcnum.Size() != surfind.Size())
                    cerr << "WARNING: solid \"" << solidname << "\" has " << surfind.Size()
                    << " surfaces and should get " << bcnum.Size() << " bc-numbers!" << endl;

                for (i = 0; i < min2(surfind.Size(),bcnum.Size()); i++)
                {
                    if (surfaces[surfind[i]]->GetBCProperty() == -1)
                        surfaces[surfind[i]]->SetBCProperty(int(bcnum[i]));
                }
            }
        }

    }

    void CSGeometry :: FindIdenticSurfaces (double eps)
    {
        int inv;
        int nsurf = GetNSurf();

        isidenticto.SetSize(nsurf);
        for (int i = 0; i < nsurf; i++)
            isidenticto[i] = i;

        for (int i = 0; i < nsurf; i++)
            for (int j = i+1; j < nsurf; j++)
            {
                if (GetSurface(j) -> IsIdentic (*GetSurface(i), inv, eps))
                {
                    INDEX_2 i2(i, j);
                    identicsurfaces.Set (i2, inv);
                    isidenticto[j] = isidenticto[i];
                }
            }
    }

    void CSGeometry ::
        GetSurfaceIndices (const Solid * sol, 
        const BoxSphere<3> & box, 
        Array<int> & locsurf) const
    {
        ReducePrimitiveIterator rpi(box);
        UnReducePrimitiveIterator urpi;

        const_cast<Solid*> (sol) -> IterateSolid (rpi);
        sol -> GetSurfaceIndices (locsurf);
        const_cast<Solid*> (sol) -> IterateSolid (urpi);

        for (int i = locsurf.Size()-1; i >= 0; i--)
        {
            bool indep = 1;
            for (int j = 0; j < i; j++)
                if (locsurf[i] == locsurf[j])
                {
                    indep = 0;
                    break;
                }

                if (!indep) locsurf.Delete(i);
        }
    }

    void CSGeometry ::
        GetIndependentSurfaceIndices (const Solid * sol, 
        const BoxSphere<3> & box, 
        Array<int> & locsurf) const
    {
        ReducePrimitiveIterator rpi(box);
        UnReducePrimitiveIterator urpi;

        ((Solid*)sol) -> IterateSolid (rpi);
        sol -> GetSurfaceIndices (locsurf);
        ((Solid*)sol) -> IterateSolid (urpi);

        for (int i = 0; i < locsurf.Size(); i++)
            locsurf[i] = isidenticto[locsurf[i]];

        for (int i = locsurf.Size()-1; i >= 0; i--)
        {
            bool indep = 1;
            for (int j = 0; j < i; j++)
                if (locsurf[i] == locsurf[j])
                {
                    indep = 0;
                    break;
                }

                if (!indep) locsurf.Delete(i);
        }

    }

    void CSGeometry ::
        GetIndependentSurfaceIndices (Array<int> & locsurf) const
    {
        for (int i = 0; i < locsurf.Size(); i++)
            locsurf[i] = isidenticto[locsurf[i]];

        for (int i = locsurf.Size()-1; i >= 0; i--)
        {
            bool indep = 1;
            for (int j = 0; j < i; j++)
                if (locsurf[i] == locsurf[j])
                {
                    indep = 0;
                    break;
                }

                if (!indep) locsurf.Delete(i);
        }
    }

    double CSGeometry ::  MaxSize () const
    {
        double maxs, mins;
        maxs = max3 (boundingbox.PMax()(0), 
            boundingbox.PMax()(1), 
            boundingbox.PMax()(2));
        mins = min3 (boundingbox.PMin()(0), 
            boundingbox.PMin()(1), 
            boundingbox.PMin()(2));
        return max2 (maxs, -mins) * 1.1;
    }

    extern CSGeometry * ParseCSG (istream & istr);

    CSGeometry*  CSGeometryRegister :: Load (string filename) const
    {
        ifstream infile(filename.c_str());

        CSGeometry *hgeom = ParseCSG (infile);
        if (!hgeom)
            throw NgException ("geo-file should start with 'algebraic3d'");

        hgeom -> FindIdenticSurfaces(1e-8f * hgeom->MaxSize()); 
        return hgeom;
    }
}
