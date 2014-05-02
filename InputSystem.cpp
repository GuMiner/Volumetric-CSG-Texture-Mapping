#include "stdafx.h"
#include "InputSystem.h"

InputSystem::ResizeEventData InputSystem::resizeEvent;

void InputSystem::Initialize(void)
{
    // Setup all the event handlers to indicate that nothing occurred.
    resizeEvent.resizeEvent = false;
}

void InputSystem::KeyTyped(GLFWwindow *pWindow, unsigned int character)
{
}

void InputSystem::KeyEvent(GLFWwindow *pWindow, int key, int scancode, int action, int mods)
{
}

void InputSystem::MouseButtonEvent(GLFWwindow *pWindow, int button, int action, int mods)
{
}

void InputSystem::ScrollEvent(GLFWwindow *pWindow, double xDelta, double yDelta)
{
}

void InputSystem::CursorTravel(GLFWwindow *pWindow, int action)
{
}

void InputSystem::CursorMove(GLFWwindow *pWindow, double xNew, double yNew)
{
}

// Simple resize handling.
void InputSystem::Resize(GLFWwindow *pWindow, int widthNew, int heightHew)
{
    resizeEvent.newWidth = widthNew;
    resizeEvent.newHeight = heightHew;
    resizeEvent.resizeEvent = true;
}
bool InputSystem::ResizeEvent(int& width, int& height)
{
    if (resizeEvent.resizeEvent)
    {
        width = resizeEvent.newWidth;
        height = resizeEvent.newHeight;
        return true;
    }

    return false;
}

// Very simple error callbacks
void InputSystem::ErrorCallback(int errCode, const char *pError)
{
    std::cout << "GLFW error " << errCode << ": " << pError << std::endl;
}
void APIENTRY InputSystem::GLCallback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar *pMessage, void *userParam)
{
    std::cout << "OpenGL debug error (" << source << ", " << type << ", " << id << ", " << severity << ", " << length << "): " << pMessage << std::endl; 
}