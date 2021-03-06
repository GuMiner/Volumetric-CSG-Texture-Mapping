Volumetric-CSG-Texture-Mapping
==============================

Synopsis
-----------------

This application is OpenGL 4.0 demonstration of [Surface-Netgen-Fork](https://github.com/GuMiner/Surface-Netgen-Fork) used for real-time CSG operations with volumetric texture mapping. Unlike a standard mesh, volumetric CSG texture msshing allows for the *holes* of an object to be moved, revealing hidden features

Status
------
This demonstration application has been completed. I have listed potential improvements at the bottom of this README file, but have no intention of performing said improvements at any time in the near future.

Overview of the Problem
-----------------------

I wanted to use constructive solid geometry (CSG) operations to build a mesh in real-time, and then display those meshes in a game. This scheme would allow for a player to shoot/cut through walls, and see/physically interact with the holes. 

For example, if the player is presented with a watermelon, the player can shoot the watermelon, revealing the inner surface and seeds. Or, the player could cut holes in the wall and later walk through them, seeing the inner structure of the wall by doing so. 

How the Problem Was Solved
--------------------------

1. Using [Surface-Netgen-Fork](https://github.com/GuMiner/Surface-Netgen-Fork), generate the mesh. Also generate a voxel-grid detailing the inner and outer structure of the object.
2. Because there are no guarantees on the triangle order returned by [Surface-Netgen-Fork](https://github.com/GuMiner/Surface-Netgen-Fork), flatten out the triangles onto a 2D texture.
3. For each triangle, determine where each pixel would lie in 3D space.
4. Because each point in 3D space maps to a voxel, color that pixel with the voxel.

Demonstration Usage Information
-------------------------------

The application opens a console window (with informational text) and a OpenGL 4.0 context. 

Key List

* 1-4: Fragment-shader mapped CSG operations of various curvature safeties.
* 5: Static File-loading example.
* 6: Triangle-identification texture mapping example.
* 7-9: Various volumetric textures mapped to the 3D surface.
* Q,A,W,S,E,D: x, y, and z camera motion.
* Spacebar: Pause (clock still runs in the background).
* Escape: Quit

Don't forget to check the console for useful information!
Also note that the executable requires the Visual Studio 2012 redistributable (x32), available from [here](http://www.microsoft.com/en-us/download/details.aspx?id=30679).

Future work
---------

0. There are noticable triangle outlines on the texture map that should be cleaned up.
1. The texture-mapping algorithm I have leaves a lot of empty space in-between triangles. Reducing this would increase the quality of the texture applied to each triangle.
2. Parts of the texture-mapping algorithm could be implemented in a separate GPU pass, speeding up the process.
3. The UI can be improved and made more user-friendly.
4. Error checking on failed calls to new -- which can happen (very rarely).
5. The triangles resulting from [Surface-Netgen-Fork](https://github.com/GuMiner/Surface-Netgen-Fork) do not have consistent windings. Either fix the generator, or rewind the triangles to allow for back-face culling.

Dependencies
-------------------
This program uses [GLEW](http://glew.sourceforge.net/) which is under the modified BSD & MIT license, [Surface-Netgen-Fork](https://github.com/GuMiner/Surface-Netgen-Fork) which is under the LGPL 2 license, and [GLFW 3](http://www.glfw.org/index.html) which is under the zlib/png license.
