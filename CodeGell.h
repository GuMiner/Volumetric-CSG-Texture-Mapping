#pragma once
#include "stdafx.h"

#include "ng\csg.hpp"

class CodeGell
{
    //// Simple vertexes ////
    struct vertex
    {
        // position
        float x;
        float y;
        float z;
        float tx;
        float ty;
    };
    struct texVertex
    {
        float tx, ty;
    };
    struct colorVertex
    {
        float r, g, b;
    };

    // General constants.
    static const int OPENGL_MAJOR = 4, OPENGL_MINOR = 0;
    static const int ALIASING_LEVEL = 1;
    static const int FPS_TARGET = 60;
    static const int TEXTURE_WH = 1024; // Texture width and height.
    const float FOV_Y, NEAR_PLANE, FAR_PLANE;

    // Color cube information and generation
    static const int GRID_SIZE = 20;
    colorVertex *pColorGrid;
    void SetupVolumetricTexture(void);

    // Window parameters
    bool running;
    bool fullscreen;
    int width, height;
    std::string title;

    // Opengl window and matrixes (CPU)
    GLFWwindow *pWindow;
    float aspect;
    gm::mat4 proj_matrix, lookAt;
    
    // Program, vertex array object, and model/projection matrixes (GPU)
    GLuint renderingProgram;
    GLuint vao;
    GLuint texture;
    GLint mv_location, proj_location, dStatus;
    
    // Veretx information
    GLuint pointBuffer;
    GLsizei vertexCount;
    
    // Volumetric texture mapping methods
    void VolumetricShadeTriangle(float *pImage, texVertex& one, texVertex& two, texVertex& three, vertex& one3, vertex& two3, vertex& three3);
    inline GLfloat IndexedPointDistance(int idxOne, int idxTwo, vertex *pPoints);
    texVertex* LayoutTriangles(vertex *pPoints, GLuint *pIndicies, int indexCount);
    void DrawLine(float* pImage, float p1x, float p1y, float p2x, float p2y);
    void ValidateBounds(gm::ivec2& vec);

    // Current camera position and look at.
    float xp, yp, zp, xk, yk, zk;

    // Variables for animation in the meshes
    double sidePos;
    bool isLeft;
    double delta;

    // Meshing variables and methods
    netgen::Mesh* GetFragmentMesh(void);
    bool meshedFromFile;
    netgen::Mesh *pFileMesh;

    // Current display mode
    enum DisplayMode {FRAGMENT_SHADING_CS8, FRAGMENT_SHADING_CS4, FRAGMENT_SHADING_CS2, FRAGMENT_SHADING_CS1, FILE_LOAD, TRIANGLE_TEXTURE_SHADING, VOLUMETRIC_TEXTURE_SHADING};
    enum VolumeMode {CHECKERBOARD, COLORGRID, WATERMELLON};
    VolumeMode volumeMode;
    DisplayMode displayMode;
    
    // Updating and rendering.
    void UpdateMesh(void);
    void UpdateOpenGL(const netgen::Mesh *pMesh);
    float GetDisplayMode(void);
    void Render(double elapsedTime);

   // Shader loading
    GLuint CompileShaders(const char rootName []);
    bool loadString(const std::string& filename, std::string& result);

    // Application setup
    void ApplicationSetup(void);
    void SetupViewport(void);
    

public:
    CodeGell(void);
    int RenderLoop(void);
    ~CodeGell(void);
};

