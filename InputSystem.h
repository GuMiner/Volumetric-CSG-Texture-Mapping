#pragma once
#include "stdafx.h"

class InputSystem
{
    typedef struct 
    {
        bool resizeEvent;
        int newWidth, newHeight;
    } ResizeEventData;
    static ResizeEventData resizeEvent;

public:
    static bool ResizeEvent(int& width, int& height);

    static void KeyTyped(GLFWwindow *pWindow, unsigned int character); // GLFWcharfun
    static void KeyEvent(GLFWwindow *pWindow, int key, int scancode, int action, int mods); // GLFWkeyfun
    static void MouseButtonEvent(GLFWwindow *pWindow, int button, int action, int mods); // GLFWmousebuttonfun
    static void ScrollEvent(GLFWwindow *pWindow, double xDelta, double yDelta); // GLFWscrollfun
    static void CursorTravel(GLFWwindow *pWindow, int action); // GLFWcursorenterfun
    static void CursorMove(GLFWwindow *pWindow, double xNew, double yNew); // GLFWcursorposfun
    static void Resize(GLFWwindow *pWindow, int widthNew, int heightHew); // GLFWwindowsizefun
    
    static void ErrorCallback(int errCode, const char *pError); //GLFWerrorfun
    static void APIENTRY GLCallback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar *pMessage, void *userParam);

    static void Initialize(void);
};

