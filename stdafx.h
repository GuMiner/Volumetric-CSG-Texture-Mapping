#pragma once

// Precompiled header file to speed up compilations

// Needed for Netgen...
#define _CRT_SECURE_NO_WARNINGS 1 

// Standard C includes
#include <cstdlib>

// Data structures
#include <vector>
#include <map>

// File and console IO
#include <iostream>
#include <fstream>

// String management
#include <sstream>
#include <string>

// Time and thread management
#include <chrono>
#include <thread>

// Miscellaneous
#include <limits>

// GLEW
#include <GL/glew.h>

// SFML
#define SFML_STATIC 1 
#include <SFML/System.hpp>
#include <SFML/Network.hpp>
#include <SFML/Audio.hpp>

// GLFW
#define GLFW_DLL 1
#include <GLFW/glfw3.h>

// General Math library.
#include "gm.h"