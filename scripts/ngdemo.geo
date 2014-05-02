# Geometry example for netgen
algebraic3d
solid cube = orthobrick(0.5, 0.5, 0.5; 1, 1, 1) and not cylinder(0.75, 0.75, -0.55; 0.75, 0.75, .1; 0.1);
solid cuber = cube and not cylinder(0.75, -0.55, 0.75; 0.75, 1.1, 0.75; 0.2);
solid sph = sphere (0.8, 0.4, 0.8; 0.2);

# Be rather careful editting this -- surfaces of geometry are very weird.
curve2d curver=(3;
    0, 0;
    0.125, 0.5;
    0.25, 0;
    2;
    2, 1, 2;
    2, 2, 3);
solid rev = revolution(0, 0.25, 0.25; 0.25, 0.25, 0.25; curver);
tlo cuber;
tlo rev;
tlo sph;