#include "stdafx.h"
#include "CodeGell.h"
#include "InputSystem.h"

// Libaries required for OpenGL, GLEW, GLFW, and SFML.
#pragma comment(lib, "opengl32")
#pragma comment(lib, "lib/glfw3dll.lib")

// All other libraries
#ifndef _DEBUG
#pragma comment(lib, "lib/glew32.lib")
#else
#pragma comment(lib, "lib/glew32d.lib")
#endif

// Does class constant and generic OpenGL setup.
CodeGell::CodeGell(void)  :
    FOV_Y (50.0f), NEAR_PLANE (0.1f), FAR_PLANE(1000.0f), pColorGrid(nullptr)
{
    running = true;
    
    // Initialize GLFW
    if (!glfwInit())
    {
        std::cout << "GLFW initialization failure!" << std::endl;
        return;
    }

    // Setup GLFW

    // Window hints
    glfwDefaultWindowHints();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, OPENGL_MAJOR);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, OPENGL_MINOR);
    glfwWindowHint(GLFW_SAMPLES, ALIASING_LEVEL);

    bool debug = false;
    if (debug)
    {
        glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);
    }
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    
    // Window size and title.
    fullscreen = false;
    width = 1366;
    height = 768;
    title = std::string("GLEW & GLFW & OpenGL Window");
    
    // Window creation.
    if (fullscreen)
    {
        pWindow  = glfwCreateWindow(width, height, title.c_str(), glfwGetPrimaryMonitor(), NULL);
    }
    else
    {
        pWindow  = glfwCreateWindow(width, height, title.c_str(), NULL, NULL);
    }

    if (pWindow == NULL)
    {
        std::cout << "GLFW window creation failure: " << glewGetErrorString(glGetError()) << std::endl;
        return;
    }
    glfwMakeContextCurrent(pWindow);
    
    // Callback initialization.
    InputSystem::Initialize();
    glfwSetCharCallback(pWindow, InputSystem::KeyTyped);
    glfwSetKeyCallback(pWindow, InputSystem::KeyEvent);
    glfwSetMouseButtonCallback(pWindow, InputSystem::MouseButtonEvent);
    glfwSetScrollCallback(pWindow, InputSystem::ScrollEvent);
    glfwSetCursorEnterCallback(pWindow, InputSystem::CursorTravel);
    glfwSetCursorPosCallback(pWindow, InputSystem::CursorMove);
    glfwSetWindowSizeCallback(pWindow, InputSystem::Resize);
    glfwSetErrorCallback(InputSystem::ErrorCallback);

    
    // GLEW initialization.
    glewExperimental = GL_TRUE;
    GLenum err = glewInit();
    if (err != GLEW_OK)
    {
        std::cout << "GLEW initialization failure: " << glewGetErrorString(err) << std::endl;
        return;
    }
    
    // Display debug info and enable debug mode if specified.
    std::cout << "Using OpenGL vendor " << glGetString(GL_VENDOR) << ", version " << glGetString(GL_VERSION) << std::endl;
    std::cout << "  >> System uses OpenGL renderer " << glGetString(GL_RENDERER) << " <<" << std::endl;
    if (debug)
    {
        if (GLEW_VERSION_4_0)
        {
            glDebugMessageCallback(InputSystem::GLCallback, nullptr);
            glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
        }
    }

    sidePos = 0.50;
    isLeft = true;
}

// Compiles and links the shaders. Based on the OpenGL SuperBible code.
GLuint CodeGell::CompileShaders(const char rootName [])
{
    GLuint program;
    GLuint vertexShader;
    GLuint fragmentShader;

    std::string vsShader, fsShader, tcsShader, tesShader, gsShader;
    std::stringstream filenameStream, filenameStream2;
    filenameStream << "shaders/" << rootName << ".vs";
    filenameStream2 << "shaders/" << rootName << ".fs";

    if (!loadString(filenameStream.str().c_str(), vsShader))
    {
        std::cout << "Could not load vertex shader!" << std::endl;
        exit(-1);
    }

    if (!loadString(filenameStream2.str().c_str(), fsShader))
    {
        std::cout << "Could not load fragment shader!" << std::endl;
        exit(-1);
    }

    const char* vss = vsShader.c_str();
    const char* fss = fsShader.c_str();

    char buffer[1024];
    GLint len;

    GLint compileStatus;
    vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vss, NULL);
    glCompileShader(vertexShader);
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &compileStatus);
    if (!compileStatus)
    {
        glGetShaderInfoLog(vertexShader, 1024, &len, buffer);
        std::cout << std::endl << "Error: " << glewGetErrorString(glGetError()) << " " << buffer << std::endl;
    }

    fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fss, NULL);
    glCompileShader(fragmentShader);
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &compileStatus);
    if (!compileStatus)
    {
        glGetShaderInfoLog(fragmentShader, 1024, &len, buffer);
        std::cout << "Error: " << glewGetErrorString(glGetError()) << " " << buffer << std::endl;
    }

    // Create program
    program = glCreateProgram();
    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);
    glLinkProgram(program);
    glGetProgramiv(program, GL_LINK_STATUS, &compileStatus);
    if (!compileStatus)
    {
        glGetProgramInfoLog(program, 1024, &len, buffer);
        std::cout << "Error: " << glewGetErrorString(glGetError()) << " " << buffer << std::endl;
    }
    
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    return program;
}

// Returns the distance between two indexed points.
inline GLfloat CodeGell::IndexedPointDistance(int idxOne, int idxTwo, vertex *pPoints)
{
    GLfloat xDistSqd = (pPoints[idxOne].x - pPoints[idxTwo].x)*(pPoints[idxOne].x - pPoints[idxTwo].x);
    GLfloat yDistSqd = (pPoints[idxOne].y - pPoints[idxTwo].y)*(pPoints[idxOne].y - pPoints[idxTwo].y);
    GLfloat zDistSqd = (pPoints[idxOne].z - pPoints[idxTwo].z)*(pPoints[idxOne].z - pPoints[idxTwo].z);
    return sqrt(xDistSqd + yDistSqd + zDistSqd);
}

// Helper method for CompileShaders
bool CodeGell::loadString(const std::string& filename, std::string& result)
{
    std::ifstream file (filename.c_str());
    if (!file) 
    {
        return false;
    }
        
    std::stringstream stream;
    stream << file.rdbuf();
    result = stream.str();
    
    file.close();
    return true;
}

// Performs a rough model unwrapping method (nearly random!) (TODO papercraft algorithms better/usable here?)
CodeGell::texVertex* CodeGell::LayoutTriangles(vertex *pPoints, GLuint *pIndicies, int indexCount)
{
    texVertex *ptexVertex = new texVertex[indexCount];

    // Figure out the bounding box of all triangles and load the texture coordinates with relative positions.
    GLfloat *pBoundingBoxes = new GLfloat [indexCount*2/3];
    GLfloat xSum = 0, ySum = 0;
    for (int i = 0; i < indexCount/3; i++)
    {
        // Triangle indicies
        GLuint idxOne   = pIndicies[i*3];
        GLuint idxTwo   = pIndicies[i*3 + 1];
        GLuint idxThree = pIndicies[i*3 + 2];

        // Distance x is the triangle max size when the one-two line is horizontal
        // Note that we don't care about the winding, as we're ok with inverted triangles for texture mapping
        GLfloat distOne = IndexedPointDistance(idxTwo, idxThree, pPoints);
        GLfloat distTwo = IndexedPointDistance(idxOne, idxThree, pPoints);
        GLfloat distThree = IndexedPointDistance(idxOne, idxTwo, pPoints);
        GLfloat oneAngle = acos((distTwo*distTwo + distThree*distThree - distOne*distOne)/(2*distTwo*distThree));

        // Both one and two y offsets at zero by definition
        ptexVertex[i*3].ty = 0;
        ptexVertex[i*3 + 1].ty = 0;

        // If the triangle is tilted to the left, apply an x-offset
        const float HALF_PI = 3.141592653589f/2.0f;
        if (oneAngle > HALF_PI)
        {

            // Same as cos(side three)*distThree, projection on the x axis.
            GLfloat cosThreeAngle = (distThree*distThree + distOne*distOne - distTwo*distTwo)/(2*distOne*distThree);

            // X offset, bottom straight line
            ptexVertex[i*3].tx = cosThreeAngle*distOne - distThree;
            ptexVertex[i*3 + 1].tx = cosThreeAngle*distOne; // Second point
            
            // By definition at zero from x-offset
            ptexVertex[i*3 + 2].tx = 0;
            ptexVertex[i*3 + 2].ty = sqrt(distOne*distOne - cosThreeAngle*cosThreeAngle*distOne*distOne); // y-projection

            // X and Y bounding boxes
            pBoundingBoxes[i*2] = cosThreeAngle*distOne;
            pBoundingBoxes[i*2 + 1] = ptexVertex[(i*3 + 2)].ty; // Maximum height from y-projection.
        }
        else // Triangle tilted to the right
        {
            pBoundingBoxes[i*2] = max(distThree, cos(oneAngle)*distTwo);
            pBoundingBoxes[i*2 + 1] = sqrt(distTwo*distTwo - cos(oneAngle)*distTwo*cos(oneAngle)*distTwo);
            
            // x positions, rather straightforward
            ptexVertex[i*3].tx = 0;
            ptexVertex[i*3 + 1].tx = distThree;
            ptexVertex[i*3 + 2].tx = cos(oneAngle)*distTwo;

            ptexVertex[i*3 + 2].ty = pBoundingBoxes[i*2 + 1];
        }

        xSum += pBoundingBoxes[i*2];
        ySum += pBoundingBoxes[i*2 + 1];
    }

    // First pass done. Now try to place and squarify the triangles by skewing the distribution.
    // std::cout << "x sum: " << xSum << ", y sum: " << ySum << std::endl;
    int trianglesPerYHeight = (int)(sqrt(indexCount/3)*xSum/ySum);
    // std::cout << "Num triangles: " << trianglesPerYHeight << std::endl;

    GLfloat xSideMax = 0;
    GLfloat yMax = 0;

    GLfloat xCurrent = 0, yCurrent = 0;
    int tCount = 0;
    for (int i = 0; i < indexCount/3; i++)
    {
        
        // Add x position offsets
        ptexVertex[i*3].tx += xCurrent;
        ptexVertex[i*3 + 1].tx += xCurrent;
        ptexVertex[i*3 + 2].tx += xCurrent;
        
        xSideMax = max(xSideMax, pBoundingBoxes[i*2]);

        // Add y position offsets
        ptexVertex[i*3].ty += yCurrent;
        ptexVertex[i*3 + 1].ty += yCurrent;
        ptexVertex[i*3 + 2].ty += yCurrent;
        yCurrent += pBoundingBoxes[i*2 + 1];

        // Update maximum y position.
        yMax = max(yMax, yCurrent);

        // Jump to the next series of triangles
        ++tCount;
        if (tCount > trianglesPerYHeight)
        {
            tCount = 0;

            xCurrent += xSideMax;
            xSideMax = 0;

            yCurrent = 0;
        }
    }

    delete [] pBoundingBoxes;

    // Update xCurrent to be xMax.
    xCurrent += xSideMax;
    
    //std::cout << "Actual maxes: (" << xCurrent << ", " << yMax << ")" << std::endl;

    // Finally, normalize placed triangles
    for (int i = 0; i < indexCount; i++)
    {
        ptexVertex[i].tx /= xCurrent;
        ptexVertex[i].ty /= yMax;
    }

    return ptexVertex;
}

// Validate that we won't exceed our array bounds.
void CodeGell::ValidateBounds(gm::ivec2& vec)
{
    if(vec[0] < 0)
    {
        vec[0] = 0;
    }
    if (vec[1] < 0)
    {
        vec[1] = 0;
    }
    if (vec[0] >= TEXTURE_WH)
    {
        vec[0] = TEXTURE_WH - 1;
    }
    if (vec[1] >= TEXTURE_WH)
    {
        vec[1] = TEXTURE_WH - 1;
    }
}

// Shades the triangle volumetrically using the color voxels.
void CodeGell::VolumetricShadeTriangle(float *pImage, texVertex& one, texVertex& two, texVertex& three, vertex& one3, vertex& two3, vertex& three3)
{
    // First, convert to integer vertices & check bounds.
    // Overall, I should really use the gm.h header and try to fixup some of the NetGen data accesses as well. 
    gm::ivec2 oneI = gm::ivec2((int)(one.tx*TEXTURE_WH + 0.5f), (int)(one.ty*TEXTURE_WH));
    gm::ivec2 twoI = gm::ivec2((int)(two.tx*TEXTURE_WH + 0.5f), (int)(two.ty*TEXTURE_WH));
    gm::ivec2 threeI = gm::ivec2((int)(three.tx*TEXTURE_WH + 0.5f), (int)(three.ty*TEXTURE_WH + 0.5f));
    ValidateBounds(oneI);
    ValidateBounds(twoI);
    ValidateBounds(threeI);
    
    float xPos = (float)oneI[0];
    float xLength = (float)(twoI[0] - oneI[0]);
    gm::vec2 yvec = gm::vec2((float)(threeI[0] - oneI[0]), (float)(threeI[1] - oneI[1])).normalize();
    gm::vec2 yvecinv = gm::vec2((float)(threeI[0] - twoI[0]), (float)(threeI[1] - twoI[1])).normalize();

    // NOT normalized.
    gm::vec3 x3v = gm::vec3(two3.x - one3.x, two3.y - one3.y, two3.z - one3.z);
    gm::vec3 y3v = gm::vec3(three3.x -  one3.x, three3.y - one3.y, three3.z - one3.z);
    
    // Go through all the layers of the triangle.
    for (int j = oneI[1]; j < threeI[1]; j++)
    {
        for (int i = (int)xPos - 1; i < (int)(xPos + xLength + 1); i++)
        {
            if (i < 0)
            {
                continue;
            }
            if (i >= TEXTURE_WH)
            {
                continue;
            }

            // Weight current position into 3D space.
            float fracy = (float)(j - oneI[1])/(float)(threeI[1] - oneI[1]);
            float fracx = ((float)(i - oneI[0]) - fracy*(float)(threeI[0] - oneI[0]))/(float)(twoI[0] - oneI[0]);
            
            // Gridify the position in 3d
            gm::vec3 pos3d = gm::vec3(x3v[0]*fracx, x3v[1]*fracx, x3v[2]*fracx) + gm::vec3(y3v[0]*fracy, y3v[1]*fracy, y3v[2]*fracy) + gm::vec3(one3.x, one3.y, one3.z);
            gm::ivec3 grids = gm::ivec3((int)((float)GRID_SIZE*pos3d[0]), (int)((float)GRID_SIZE*pos3d[1]), (int)((float)GRID_SIZE*pos3d[2]));
            if (grids[0] < 0){ grids[0] = 0; } if (grids[1] < 0) { grids[1] = 0; } if (grids[2] < 0) { grids[2] = 0; }
            if (grids[0] >= GRID_SIZE) { grids[0] = GRID_SIZE - 1; } if (grids[1] >= GRID_SIZE) { grids[1] = GRID_SIZE - 1; } if (grids[2] >= GRID_SIZE) { grids[2] = GRID_SIZE - 1; }

            // Assign the volumetric color (non-aliased)
            pImage[(j*TEXTURE_WH + i)*4] = pColorGrid[grids[2]*GRID_SIZE*GRID_SIZE + grids[1]*GRID_SIZE + grids[0]].r;
            pImage[(j*TEXTURE_WH + i)*4 + 1] = pColorGrid[grids[2]*GRID_SIZE*GRID_SIZE + grids[1]*GRID_SIZE + grids[0]].g;
            pImage[(j*TEXTURE_WH + i)*4 + 2] = pColorGrid[grids[2]*GRID_SIZE*GRID_SIZE + grids[1]*GRID_SIZE + grids[0]].b;
            pImage[(j*TEXTURE_WH + i)*4 + 3] = 1.0f;
        }
        xPos += (yvec[0]/yvec[1]);
        xLength += ((yvecinv[0]/yvecinv[1]) - (yvec[0]/yvec[1]));
    }
}

// Draws a line on the texture to identify triangle bounds.
void CodeGell::DrawLine(float *pImage, float p1x, float p1y, float p2x, float p2y)
{
    float dist = sqrt((p1x - p2x)*(p1x - p2x) + (p1y - p2y)*(p1y - p2y));
    int steps = (int)(dist*(float)TEXTURE_WH);
    float posx = p1x;
    float posy = p1y;
    
    float stepSizeX = (p2x - p1x)/(float)steps;
    float stepSizeY = (p2y - p1y)/(float)steps;
    for(int i = 0; i < steps; i++)
    {
        int xp = (int)(posx*TEXTURE_WH);
        int yp = (int)(posy*TEXTURE_WH);
        if (xp < 0)
        {
            xp = 0;
        }
        if (yp < 0)
        {
            yp = 0;
        }
        if (xp >= TEXTURE_WH)
        {
            xp = TEXTURE_WH - 1;
        }
        if (yp >= TEXTURE_WH)
        {
            yp = TEXTURE_WH - 1;
        }

        pImage[(yp*TEXTURE_WH + xp)*4] = 1.0f;
        pImage[(yp*TEXTURE_WH + xp)*4 + 1] = 1.0f;
        pImage[(yp*TEXTURE_WH + xp)*4 + 2] = 1.0f;
        pImage[(yp*TEXTURE_WH + xp)*4 + 3] = 1.0f;

        posx += stepSizeX;
        posy += stepSizeY;
    }
}

// Updates the OpenGL display from the contents of the mesh.
void CodeGell::UpdateOpenGL(const netgen::Mesh *pMesh)
{
     /**** Note that we need to dynamically-allocate pretty much everything, because of the the mesh changes! *****/
    texVertex *pTexVertex = nullptr;

    // Copy in the points
    int pointCount = pMesh->GetNP()*3;
    vertex *pVertices = new vertex [pointCount];
    for (int i = 1; i <= pMesh->GetNP(); i++)
    {
        const netgen::Point3d& p = pMesh->Point(i);
        pVertices[i - 1].x = (float)p.X();
        pVertices[i - 1].y = (float)p.Y();
        pVertices[i - 1].z = (float)p.Z();
    }
    
    // Copy in the indicies
    int indexCount = pMesh->GetNSE()*3;
    GLuint *pIndicies = new GLuint[indexCount];
    for (int i = 1; i <= pMesh->GetNSE(); i++)
    {
        const netgen::Element2d& el = pMesh->SurfaceElement(i);

        for (int j = 1; j <= el.GetNP(); j++)
        {
            pIndicies[(i - 1)*3 + (j - 1)] = (el.PNum(j) - 1);
        }
    }

    /** We don't delete yet, as we now need to determine UV mapping and ID our triangles! **/

    // Adds the triangle ID or volumetric texture mapping if desired.
    if (displayMode == DisplayMode::TRIANGLE_TEXTURE_SHADING || displayMode == DisplayMode::VOLUMETRIC_TEXTURE_SHADING)
    {
        // First, layout our triangles. These vertices are ordered for each vertex.
        pTexVertex = LayoutTriangles(pVertices, pIndicies, indexCount);

        // Now draw dots at our texture coordinates (for now)
        int idSize = TEXTURE_WH*TEXTURE_WH*4;
        float *pImageData = new float[idSize];
        for (int i = 0; i < idSize/4; i++)
        {
            pImageData[i*4] = 0.25f;
            pImageData[i*4 + 1] = 0.25f;
            pImageData[i*4 + 2] = 0.25f;
            pImageData[i*4 + 3] = 1.0f;
        }


        // Draw colored lines for all the triangles to verify their existence.
        for (int i = 0; i < indexCount/3; i++)
        {
            if (displayMode == DisplayMode::TRIANGLE_TEXTURE_SHADING)
            {
                DrawLine(pImageData, pTexVertex[i*3].tx, pTexVertex[i*3].ty, pTexVertex[i*3 + 1].tx, pTexVertex[i*3 + 1].ty);
                DrawLine(pImageData, pTexVertex[i*3 + 2].tx, pTexVertex[i*3 + 2].ty, pTexVertex[i*3 + 1].tx, pTexVertex[i*3 + 1].ty);
                DrawLine(pImageData, pTexVertex[i*3].tx, pTexVertex[i*3].ty, pTexVertex[i*3 + 2].tx, pTexVertex[i*3 + 2].ty);
            }
            else if (displayMode == DisplayMode::VOLUMETRIC_TEXTURE_SHADING)
            {
                VolumetricShadeTriangle(pImageData, pTexVertex[i*3], pTexVertex[i*3 + 1], pTexVertex[i*3 + 2], pVertices[pIndicies[i*3]], pVertices[pIndicies[i*3 + 1]], pVertices[pIndicies[i*3 + 2]]);
            }
        }

        // Pass in the texture to OpenGL
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, TEXTURE_WH, TEXTURE_WH, GL_RGBA, GL_FLOAT, pImageData);

        delete [] pImageData;
    }

    // Interleave our data (and convert to non-indexed form) and pass it into OpenGL.
    vertex *pGLVertices = new vertex[indexCount];
    for (int i = 0; i < indexCount; i++)
    {
        // Lookup the vertex coordinates.
        int vertex = pIndicies[i];
        pGLVertices[i].x = pVertices[vertex].x;
        pGLVertices[i].y = pVertices[vertex].y;
        pGLVertices[i].z = pVertices[vertex].z;
        
        // Texture coordinates will already be in sync from lookup
        if (pTexVertex != nullptr)
        {
            pGLVertices[i].tx = pTexVertex[i].tx;
            pGLVertices[i].ty = pTexVertex[i].ty;
        }
    }

    // Update the OpenGL buffers.
    glBufferData(GL_ARRAY_BUFFER, indexCount*sizeof(vertex), pGLVertices, GL_DYNAMIC_DRAW);
    vertexCount = indexCount; // Well, now it does...

    // Done with the temporary indicies and points.
    delete [] pIndicies;
    delete [] pVertices;
    delete [] pGLVertices;

    if (pTexVertex != nullptr)
    {
        delete [] pTexVertex;
    }
}

// Generates a cube with two cylinders cut out of it that is animated (based on the refresh rate).
netgen::Mesh* CodeGell::GetFragmentMesh(void)
{
    // Update the position of the movable pieces for the fragment mesh.
    sidePos += isLeft ? delta : -delta;
    if (sidePos > 0.60f)
    {
        isLeft = false;
        sidePos -= delta;
    }
    else if (sidePos < 0.40f)
    {
        isLeft = true;
        sidePos += delta;
    }

    // Generate the realtime fragment mesh.
    netgen::CSGeometry *pGeom = new netgen::CSGeometry;
    
    netgen::Point<3> pa = netgen::Point<3>(0, 0, 0);
    netgen::Point<3> pb = netgen::Point<3>(1, 1, 1);
    netgen::Primitive *ncube = new netgen::OrthoBrick(pa, pb);
    pGeom->AddSurfaces(ncube);
    netgen::Solid *cube = new netgen::Solid(ncube);

    netgen::Point<3> pac = netgen::Point<3>(0.5, sidePos, -1.1);
    netgen::Point<3> pbc = netgen::Point<3>(0.5, sidePos, 1.1);
    double r = 0.3;
    netgen::OneSurfacePrimitive *surf = new netgen::Cylinder(pac, pbc, r);
    pGeom->AddSurfaces(surf);
    netgen::Solid *cyl = new netgen::Solid(surf);

    netgen::Point<3> pad = netgen::Point<3>(0.5, -1.1, sidePos);
    netgen::Point<3> pbd = netgen::Point<3>(0.5, 1.1, sidePos);
    double r2 = 0.2;
    netgen::OneSurfacePrimitive *surfr = new netgen::Cylinder(pad, pbd, r2);
    pGeom->AddSurfaces(surfr);
    netgen::Solid *cylr = new netgen::Solid(surfr);

    netgen::Solid *cylinder = new netgen::Solid(netgen::Solid::SUB, cyl);
    netgen::Solid  *solid = new netgen::Solid(netgen::Solid::SECTION, new netgen::Solid(netgen::Solid::SECTION, cube, cylinder), new netgen::Solid(netgen::Solid::SUB, cylr));
    netgen::Flags flags;

    pGeom->SetSolid("Root", new netgen::Solid(netgen::Solid::ROOT, solid));
    pGeom->SetTopLevelObject(solid);
    pGeom->SetFlags("Root", flags);
    pGeom->FindIdenticSurfaces(1e-8*pGeom->MaxSize()); // Initializes surfaces

    netgen::Mesh *pMesh = NULL;
    netgen::MeshingParameters mParams;
    switch(displayMode)
    {
    case DisplayMode::FRAGMENT_SHADING_CS1:
        mParams.curvaturesafety = 1;
        break;
    case DisplayMode::FRAGMENT_SHADING_CS2:
        mParams.curvaturesafety = 2;
        break;
    case DisplayMode::FRAGMENT_SHADING_CS4:
        mParams.curvaturesafety = 4;
        break;
    case DisplayMode::FRAGMENT_SHADING_CS8:
        mParams.curvaturesafety = 8;
        break;
    default: break; // Default of 8 curvature safety.
    }
    pGeom->GenerateMesh(pMesh, mParams); // Runs from analyze to surface meshing automatically.
    
    // Done with the geometry.
    delete pGeom;

    return pMesh;
}

// Updates the meshes used in the application.
void CodeGell::UpdateMesh(void)
{
    // Does a static file load for high-speed non-changing objects.
    if (displayMode == DisplayMode::FILE_LOAD)
    {
        // Run the file mesh if we haven't already done so.
        if (!meshedFromFile)
        {
            std::string filename ("scripts/ngdemo.geo");
            netgen::CSGeometryRegister regr;
            netgen::CSGeometry *pGeom = regr.Load(filename);

            netgen::MeshingParameters mParams;
            pGeom->GenerateMesh(pFileMesh, mParams);

            delete pGeom;
            meshedFromFile = true;
        }
        
        UpdateOpenGL(pFileMesh);
        return;
    }
    
    // Gets the dynamically-generated CSG mesh.
    netgen::Mesh *pMesh = GetFragmentMesh();
    if (pMesh == NULL)
    {
        std::cout << "Something really bad happened. Exiting" << std::endl;
        exit(-1);
    }

    UpdateOpenGL(pMesh);

    // Done with the mesh.
    delete pMesh;
}

// Sets up a grid that acts as a volumetric texture, mapped as a color box.
void CodeGell::SetupVolumetricTexture(void)
{
    // May occur when switching modes.
    if (pColorGrid != nullptr)
    {
        delete [] pColorGrid;
    }

    pColorGrid = new colorVertex[GRID_SIZE*GRID_SIZE*GRID_SIZE];
    for (int i = 0; i < GRID_SIZE; i++)
    {
        for (int j = 0; j < GRID_SIZE; j++)
        {
            for (int k = 0; k < GRID_SIZE; k++)
            {
                if (volumeMode == VolumeMode::CHECKERBOARD)
                {
                    // B&W grid
                    pColorGrid[k*GRID_SIZE*GRID_SIZE + j*GRID_SIZE + i].r = (float)((i + j + k) % 2);
                    pColorGrid[k*GRID_SIZE*GRID_SIZE + j*GRID_SIZE + i].g = (float)((i + j + k) % 2);
                    pColorGrid[k*GRID_SIZE*GRID_SIZE + j*GRID_SIZE + i].b = (float)((i + j + k) % 2);
                }
                else if (volumeMode == VolumeMode::COLORGRID)
                {
                    // Color grid
                    pColorGrid[k*GRID_SIZE*GRID_SIZE + j*GRID_SIZE + i].r = (float)i/(float)GRID_SIZE;
                    pColorGrid[k*GRID_SIZE*GRID_SIZE + j*GRID_SIZE + i].g = (float)j/(float)GRID_SIZE;
                    pColorGrid[k*GRID_SIZE*GRID_SIZE + j*GRID_SIZE + i].b = (float)k/(float)GRID_SIZE;
                }
                else if (volumeMode == VolumeMode::WATERMELLON)
                {
                    // Watermelon
                    if (i == 0 || j == 0 || k == 0 || i == GRID_SIZE - 1 || j == GRID_SIZE -1 || k == GRID_SIZE - 1)
                    {
                        // Skin
                        pColorGrid[k*GRID_SIZE*GRID_SIZE + j*GRID_SIZE + i].r = 34.f/255.0f;
                        pColorGrid[k*GRID_SIZE*GRID_SIZE + j*GRID_SIZE + i].g = 177.0f/255.0f;
                        pColorGrid[k*GRID_SIZE*GRID_SIZE + j*GRID_SIZE + i].b = 76.0f/255.0f;
                    }
                    else
                    {
                        // Seeds and inner
                        bool isSeed = rand() % 4 == 0;
                        
                        if (isSeed)
                        {
                            float color = 82.0f/255.0f;
                            pColorGrid[k*GRID_SIZE*GRID_SIZE + j*GRID_SIZE + i].r = color;
                            pColorGrid[k*GRID_SIZE*GRID_SIZE + j*GRID_SIZE + i].g = color;
                            pColorGrid[k*GRID_SIZE*GRID_SIZE + j*GRID_SIZE + i].b = color;
                        }
                        else
                        {
                            pColorGrid[k*GRID_SIZE*GRID_SIZE + j*GRID_SIZE + i].r = 1.0f;
                            pColorGrid[k*GRID_SIZE*GRID_SIZE + j*GRID_SIZE + i].g = 0.5f;
                            pColorGrid[k*GRID_SIZE*GRID_SIZE + j*GRID_SIZE + i].b = 0.5f;
                        }
                    }
                }
            }
        }
    }
}

// Perform once-through initialization, such as enabling and generating the buffers, viewport, depth test, etc...
void CodeGell::ApplicationSetup(void)
{
    // OpenGL setup
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // Setup of vertex transfer (note we're using the "vertex" object in CodeGell)
    glGenBuffers(1, &pointBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, pointBuffer);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), (GLvoid*)offsetof(vertex, x));
    glEnableVertexAttribArray(0);

    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(vertex), (GLvoid*)offsetof(vertex, tx));
    glEnableVertexAttribArray(1);

    // Setup texture and tc coordinate transfer
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA32F, TEXTURE_WH, TEXTURE_WH);

    // Setup of shaders
    renderingProgram = CompileShaders("render");
    mv_location = glGetUniformLocation(renderingProgram, "mv_matrix");
    proj_location = glGetUniformLocation(renderingProgram, "proj_matrix");
    dStatus = glGetUniformLocation(renderingProgram, "dStatus");
    
    // Only works if faces are positioned appropriately, which NetGen does NOT do. Another TODO...
    // glEnable(GL_CULL_FACE);
    // glFrontFace(GL_CW);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    
    SetupViewport();

    // Setup view position
    xp = 0; 
    yp = 0;
    zp = 0;
    xk = xp;
    yk = yp;
    zk = zp + 6;

    // Setup animation and item display parameters.
    meshedFromFile = false;
    pFileMesh = NULL;
    delta = 0.05f;
    displayMode = DisplayMode::FRAGMENT_SHADING_CS2;
    volumeMode = VolumeMode::CHECKERBOARD;

    SetupVolumetricTexture();
}

// Perform shutdown tasks.
CodeGell::~CodeGell(void)
{
    // Destroy mesh if from file.
    if (pFileMesh != NULL)
    {
        delete pFileMesh;
    }

    delete [] pColorGrid;

    // Application shutdown
    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(1, &pointBuffer);
    glDeleteTextures(1, &texture);

    glDeleteProgram(renderingProgram);

    // Windowing shutdown
    glfwDestroyWindow(pWindow);
    glfwTerminate();
}

// Sizes the viewport as appropriate.
void CodeGell::SetupViewport(void)
{
    // Viewing aspect ratio and projection matrix.
    aspect = (float)width/ (float)height;
    proj_matrix = gm::perspective(FOV_Y, aspect, NEAR_PLANE, FAR_PLANE);
    glViewport(0, 0, width, height);
}

// Run the application.
int CodeGell::RenderLoop(void)
{
    ApplicationSetup();

    double timeDelta = 1.0f/(double)FPS_TARGET;
    double lastTime = (double)glfwGetTime();
    while (running)
    {
        // Update system
        lookAt = gm::lookat(gm::vec3(xp, yp, zp), gm::vec3(xp, yp, zp + 6), gm::vec3(0, 1, 0));

        glUseProgram(renderingProgram);

        UpdateMesh();

        // Draw and swap buffers
        Render((double)glfwGetTime());
        glfwSwapBuffers(pWindow);
       
        while (glfwGetKey(pWindow, GLFW_KEY_SPACE))
        {
            glfwPollEvents();
        }

        // Handle events.
        glfwPollEvents();
        if (glfwGetKey(pWindow, GLFW_KEY_ESCAPE) == GLFW_PRESS || glfwWindowShouldClose(pWindow))
        {
            running = false;
        }

        // Handle window resizes
        if (InputSystem::ResizeEvent(width, height))
        {
            SetupViewport();
        }

        // Handle key presses for motion
        bool motionPrintout = false;
        if (glfwGetKey(pWindow, GLFW_KEY_Q))
        {
            xp += 0.1f;
            motionPrintout = true;
        }
        if (glfwGetKey(pWindow, GLFW_KEY_A))
        {
            xp -= 0.1f;
            motionPrintout = true;
        }
        if (glfwGetKey(pWindow, GLFW_KEY_W))
        {
            yp += 0.1f;
            motionPrintout = true;
        }
        if (glfwGetKey(pWindow, GLFW_KEY_S))
        {
            yp -= 0.1f;
            motionPrintout = true;
        }
        if (glfwGetKey(pWindow, GLFW_KEY_E))
        {
            zp += 0.1f;
            motionPrintout = true;
        }
        if (glfwGetKey(pWindow, GLFW_KEY_D))
        {
            zp -= 0.1f;
            motionPrintout = true;
        }

        if (glfwGetKey(pWindow, GLFW_KEY_T))
        {
            xk += 0.1f;
            motionPrintout = true;
        }
        if (glfwGetKey(pWindow, GLFW_KEY_G))
        {
            xk -= 0.1f;
            motionPrintout = true;
        }
        if (glfwGetKey(pWindow, GLFW_KEY_Y))
        {
            yk += 0.1f;
            motionPrintout = true;
        }
        if (glfwGetKey(pWindow, GLFW_KEY_H))
        {
            yk -= 0.1f;
            motionPrintout = true;
        }
        if (glfwGetKey(pWindow, GLFW_KEY_U))
        {
            zk += 0.1f;
            motionPrintout = true;
        }
        if (glfwGetKey(pWindow, GLFW_KEY_J))
        {
            zk -= 0.1f;
            motionPrintout = true;
        }

        if (motionPrintout)
        {
            std::cout << "At (" << xp << ", " << yp << ", " << zp << ") looking at (" << xk << ", " << yk << ", " << zk << ")." << std::endl;
        }

        // Handle key presses for mode changes.
        if (glfwGetKey(pWindow, GLFW_KEY_1))
        {
            displayMode = DisplayMode::FRAGMENT_SHADING_CS8;
            std::cout << "Switched to fragment shading w/ curvature safety of 8." << std::endl;
        }
        else if (glfwGetKey(pWindow, GLFW_KEY_2))
        {
            displayMode = DisplayMode::FRAGMENT_SHADING_CS4;
            std::cout << "Switched to fragment shading w/ curvature safety of 4." << std::endl;
        }
        else if (glfwGetKey(pWindow, GLFW_KEY_3))
        {
            displayMode = DisplayMode::FRAGMENT_SHADING_CS2;
            std::cout << "Switched to fragment shading w/ curvature safety of 2." << std::endl;
        }
        else if (glfwGetKey(pWindow, GLFW_KEY_4))
        {
            displayMode = DisplayMode::FRAGMENT_SHADING_CS1;
            std::cout << "Switched to fragment shading w/ curvature safety of 1." << std::endl;
        }
        else if (glfwGetKey(pWindow, GLFW_KEY_5))
        {
            displayMode = DisplayMode::FILE_LOAD;
            std::cout << "Switched to static file-load method." << std::endl;
        }
        else if (glfwGetKey(pWindow, GLFW_KEY_6))
        {
            displayMode = DisplayMode::TRIANGLE_TEXTURE_SHADING;
            std::cout << "Switched to triangle ID texture shading." << std::endl;
        }
        else if (glfwGetKey(pWindow, GLFW_KEY_7))
        {
            displayMode = DisplayMode::VOLUMETRIC_TEXTURE_SHADING;
            volumeMode = VolumeMode::CHECKERBOARD;
            std::cout << "Switched to gradient-based volumentric texture shading, checkerboard." << std::endl;;
            SetupVolumetricTexture();
        }
        else if (glfwGetKey(pWindow, GLFW_KEY_8))
        {
            displayMode = DisplayMode::VOLUMETRIC_TEXTURE_SHADING;
            volumeMode = VolumeMode::COLORGRID;
            std::cout << "Switched to gradient-based volumentric texture shading, colorgrid." << std::endl;;
            SetupVolumetricTexture();
        }
        else if (glfwGetKey(pWindow, GLFW_KEY_9))
        {
            displayMode = DisplayMode::VOLUMETRIC_TEXTURE_SHADING;
            volumeMode = VolumeMode::WATERMELLON;
            std::cout << "Switched to gradient-based volumentric texture shading, watermellon." << std::endl;;
            SetupVolumetricTexture();
        }

        // Update timer and try to sleep for the FPS Target.
        timeDelta = (double)glfwGetTime() - lastTime;
        lastTime  = (double)glfwGetTime();

        std::chrono::milliseconds sleepTime ((int)(1.0/(double)FPS_TARGET - 1000*timeDelta));
        if (sleepTime > std::chrono::milliseconds(0))
        {
            std::this_thread::sleep_for(sleepTime);
        }
    }

    return 0;
}

// Returns the current display mode
float CodeGell::GetDisplayMode(void)
{
    if (displayMode == DisplayMode::TRIANGLE_TEXTURE_SHADING || displayMode == DisplayMode::VOLUMETRIC_TEXTURE_SHADING)
    {
        return 1.0f;
    }

    return 0.0f;
}

// Draws the application!
void CodeGell::Render(double curTime)
{ 
    const GLfloat  color[] = {0, 0, 0, 1};
    const GLfloat  one = 1.0f;
    glClearBufferfv(GL_COLOR, 0, color);
    glClearBufferfv(GL_DEPTH, 0, &one);

    glUniformMatrix4fv(proj_location, 1, GL_FALSE, proj_matrix*lookAt);
    glUniform1f(dStatus, GetDisplayMode());

    // Draw the spinning cubes
    //glPointSize(10);
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            for (int k = 0; k < 1; k++)
            {
                gm::mat4 mv_matrix = gm::rotate(0.0f, 0.0f, 1.0f, 0.0f);//gm::rotate(5*(float)curTime, 0.0f, 1.0f, 0.0f);
                glUniformMatrix4fv(mv_location, 1, GL_FALSE, mv_matrix);
                glDrawArraysInstanced(GL_TRIANGLES, 0, vertexCount, 100);
            }
        }
    }
}

// Runs the main application.
int main(int argc, char* argv [])
{
    CodeGell *pGell = new CodeGell();
    int runStatus = pGell->RenderLoop();
    delete pGell;

    // Wait before closing for display purposes.
    std::cout << std::endl << "End of application. " << std::endl;
    std::chrono::milliseconds sleepTime(1000);
    std::this_thread::sleep_for(sleepTime);

    return runStatus;
}