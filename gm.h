#pragma once

#include <cstdlib>
#include <cmath>

// General Math namespace
//  Contains float ing-point specific math functions.
//  Templates are only used to avoid repetition -- else float ing-point values are preferred.
//  Naming follows GLSL specifications for vector and matrix names.
namespace gm
{
    // Angle conversion routines
    const float  PI = 3.141592635389f;

    inline float  ToDegrees(float  radAngle)
    {
        return radAngle*180.0f/PI;
    }

    inline float  ToRadians(float  degAngle)
    {
        return degAngle*PI/180.0f;
    }

    // float ing-point random numbers.
    inline float  Random()
    {
        return (float )rand()/(float )RAND_MAX;
    }

    // Vectors
    // Including the typename templage to allow for integer and float ing point vectors.
    template <typename T, const int length> 
    class vecX
    {
    protected:

        // Vector position data.
        T data[length];
    public:

        // Data element accessing operators
        inline T& operator [] (int n) 
        {
            return data[n];
        }
        inline const T& operator [] (int n) const 
        {
            return data[n];
        }

        // Default, copy, and assignment constructors
        inline vecX()
        { }
        inline vecX(const vecX& other)
        {
            for (int n = 0; n < length; n++)
            {
                data[n] = other[n];
            }
        }
        inline vecX& operator=(const vecX& other)
        {
            for (int n = 0; n < length; n++)
            {
                data[n] = other[n];
            }
            return *this;
        }

        // Construction from a scalar
        inline vecX(T other)
        {
            for (int n = 0; n < length; n++)
            {
                data[n] = other;
            }
        }

        // Dereferencing operator
        inline operator const T* () const
        {
            return &data[0];
        }

        inline static int size()
        {
            return length;
        }

        // Basic mathematical operators

        // Dot product
        inline T dot (const vecX<T, length>& other) const
        {
            T result;
            for (int n = 0; n < length; n++)
            {
                result += (data[n]*other[n]);
            }

            return result;
        }

        // Vector length
        inline T len() const
        {
            T totLen = 0;
            for (int n = 0; n < length; n++)
            {
                totLen += data[n]*data[n];
            }

            return static_cast<T>(sqrt(totLen));
        }

        // Normalization
        inline vecX normalize() const
        {
            class vecX<T, length> result;
            T totLen = len();

            for (int n = 0; n < length; n++)
            {
                result[n] = data[n]/totLen;
            }

            return result;
        }

        // Distance between two vectors
        inline T distance(const vecX<T, length>& other)
        {
            return (*this - other).len();
        }

        // Angle between two vectors
        inline T angle(const vecX<T, length>& other)
        {
            return static_cast<T>(acos(dot(other)));
        }

        // Memberwise division.
        inline vecX& operator / (const vecX& other) const
        {
            class vecX<T, length> result;
            for (int n = 0; n < length; n++)
            {
                result[n] = data[n] / other[n];
            }

            return result;
        }
        inline vecX& operator /= (const vecX& other)
        {
            for (int n = 0; n < length; n++)
            {
                data[n] /= other[n];
            }

            return *this;
        }

        // Single-element division
        inline vecX operator / (const T& other) const
        {
            class vecX<T, length> result;
            for (int n = 0; n < length; n++)
            {
                result[n] = data[n] / other;
            }

            return result;
        }
        inline vecX& operator /= (const T& other)
        {
            for (int n = 0; n < length; n++)
            {
                data[n] /= other;
            }

            return *this;
        }

        // Memberwise multiplication
        inline vecX& operator * (const vecX& other) const
        {
            class vecX<T, length> result;
            for (int n = 0; n < length; n++)
            {
                result[n] = data[n] * other[n];
            }

            return result;
        }
        inline vecX& operator *= (const vecX& other)
        {
            for (int n = 0; n < length; n++)
            {
                data[n] *= other[n];
            }

            return *this;
        }

        // Single-element multiplication
        inline vecX operator * (const T& other) const
        {
            class vecX<T, length> result;
            for (int n = 0; n < length; n++)
            {
                result[n] = data[n] * other;
            }

            return result;
        }
        inline vecX& operator *= (const T& other)
        {
            for (int n = 0; n < length; n++)
            {
                data[n] *= other;
            }

            return *this;
        }

        // Memberwise addition
        inline vecX operator + (const vecX& other) const
        {
            class vecX<T, length> result;
            for (int n = 0; n < length; n++)
            {
                result[n] = data[n] + other[n];
            }

            return result;
        }
        inline vecX& operator += (const vecX& other)
        {
            for (int n = 0; n < length; n++)
            {
                data[n] += other[n];
            }

            return *this;
        }

        // Single-element addition
        inline vecX operator + (const T& other) const
        {
            class vecX<T, length> result;
            for (int n = 0; n < length; n++)
            {
                result[n] = data[n] + other;
            }

            return result;
        }
        inline vecX& operator += (const T& other)
        {
            for (int n = 0; n < length; n++)
            {
                data[n] += other;
            }

            return *this;
        }

        // Memberwise subtraction
        inline vecX operator - (const vecX& other) const
        {
            class vecX<T, length> result;
            for (int n = 0; n < length; n++)
            {
                result[n] = data[n] - other[n];
            }

            return result;
        }
        inline vecX& operator -= (const vecX& other)
        {
            for (int n = 0; n < length; n++)
            {
                data[n] -= other[n];
            }

            return *this;
        }

        // Single-element subtraction
        inline vecX operator - (const T& other) const
        {
            class vecX<T, length> result;
            for (int n = 0; n < length; n++)
            {
                result[n] = data[n] - other;
            }

            return result;
        }
        inline vecX& operator -= (const T& other)
        {
            for (int n = 0; n < length; n++)
            {
                data[n] -= other;
            }

            return *this;
        }
        
        // Vector negation
        inline vecX operator-() const
        {
            class vecX<T, length> result;
            for (int n = 0; n < length; n++)
            {
                result[n] = -data[n];
            }
            return result;
        }

    };

    // Two-element vector
    template <typename T>
    class Tvec2 : public vecX<T, 2>
    {
    public:

        // Construction methods
        inline Tvec2 ()
        { }
        inline Tvec2 (const vecX<T, 2>& other) : vecX<T, 2>(other)
        { }
        inline Tvec2(T x, T y)
        {
            vecX<T, 2>::data[0] = x;
            vecX<T, 2>::data[1] = y;
        }
    };

    // Three-element vector
    template <typename T>
    class Tvec3 : public vecX<T, 3>
    {
    public:

        // Construction methods
        inline Tvec3 ()
        { }
        inline Tvec3 (const vecX<T, 3>& other) : vecX<T, 3>(other)
        { }
        inline Tvec3(T x, T y, T z)
        {
            vecX<T, 3>::data[0] = x;
            vecX<T, 3>::data[1] = y;
            vecX<T, 3>::data[2] = z;
        }

        // Cross product (defined only for 3-element vectors)
        inline Tvec3<T> cross(const vecX<T, 3>& other) const
        {
            Tvec3<T> result;
            result[0] = data[1]*other[2] - data[2]*other[1];
            result[1] = data[2]*other[0] - data[0]*other[2];
            result[2] = data[0]*other[1] - data[1]*other[0];

            return result;
        }
    };

    // Four-element vector
    template <typename T>
    class Tvec4 : public vecX<T, 4>
    {
    public:

        // Construction methods
        inline Tvec4 ()
        { }
        inline Tvec4 (const vecX<T, 4>& other) : vecX<T, 4>(other)
        { }
        inline Tvec4(T x, T y, T z, T w)
        {
            vecX<T, 4>::data[0] = x;
            vecX<T, 4>::data[1] = y;
            vecX<T, 4>::data[2] = z;
            vecX<T, 4>::data[3] = w;
        }
    };

    // float ing-point and integer vector specializations (following GLSL)
    typedef vecX<float , 1> vec1;
    typedef vecX<int, 1> ivec1;

    typedef Tvec2<float > vec2;
    typedef Tvec2<int> ivec2;

    typedef Tvec3<float > vec3;
    typedef Tvec3<int> ivec3;

    typedef Tvec4<float > vec4;
    typedef Tvec4<int> ivec4;

    // Prototype definition due to matrix use inside quaternions.
    template <typename T, const int width, const int height> class matXY;

    // Four-element quaternion
    template <typename T>
    class Tquaternion
    {
    public:
        Tvec4<T> data;

        // Constructors
        inline Tquaternion()
        { }
        inline Tquaternion(T x, T y, T z)
        {
            data[0] = x;
            data[1] = y;
            data[2] = z;
            data[3] = 0;
        }
        inline Tquaternion(T x, T y, T z, T w)
        {
            data[0] = x;
            data[1] = y;
            data[2] = z;
            data[3] = w;
        }
        inline Tquaternion(const Tquaternion& other)
        {
            data = other.data;
        }
        inline Tquaternion(const Tvec4<T>& other)
        {
            data = other;
        }

        // Data access operators
        inline T& operator [] (int n)
        {
            return data[n];
        }
        inline const T& operator [] (int n) const
        {
            return data[n];
        }

        // Basic mathematical operators
    
        // Negation
        inline Tquaternion operator - () const
        {
            return Tquaternion(-data);
        }

        // Addition
        inline Tquaternion operator + (const Tquaternion& other) const
        {
            return Tquaternion(data + other.data);
        }
        inline Tquaternion& operator += (const Tquaternion& other)
        {
            data += other.data;
            return *this;
        }

        // Subtraction
        inline Tquaternion operator - (const Tquaternion& other) const
        {
            return Tquaternion(data - other.data);
        }
        inline Tquaternion& operator -= (const Tquaternion& other)
        {
            data -= other.data;
            return *this;
        }

        // Memberwise scalar multiplication
        inline Tquaternion operator * (const T other) const
        {
            return Tquaternion(data*other);
        }
        inline Tquaternion& operator *= (const T other)
        {
            data *= other.data;
            return *this;
        }

        // Quaternion multiplication
        inline Tquaternion operator * (const Tquaternion& other) const
        {
            // Pulled straignt from vmath.h due to potential errors in the OGLSB.
            const T x1 = data[0];
            const T y1 = data[1];
            const T z1 = data[2];
            const T w1 = data[3];
            const T x2 = other[0];
            const T y2 = other[1];
            const T z2 = other[2];
            const T w2 = other[3];

            return Tquaternion(w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                               w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2,
                               w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2,
                               w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2);
        }

        // Memberwise scalar division.
        inline Tquaternion operator / (const T other) const
        {
            return Tquaternion(data/other);
        }
        inline Tquaternion& operator /= (const T other)
        {
            data /= other.data;
            return *this;
        }

        // Length calculation
        inline T len () const
        {
            return data.len();
        }

        // Normalization
        inline Tquaternion& normalize() const
        {
            return data.normalize();
        }

        // Quaternion access as a matrix
        // Pulled partially from vmath.h
        inline matXY<T, 4, 4> asMatrix() const
        {
            matXY<T,4,4> result;
            
            result[0][0] = T(1) - T(2) * (y*y + z*z);
            result[0][1] = T(2) * (x*y - z*w);
            result[0][2] = T(2) * (x*z + y*w);
            result[0][3] = T(0);

            result[1][0] = T(2) * (xy + z*w);
            result[1][1] = T(1) - T(2) * (x*x + z*z);
            result[1][2] = T(2) * (y*z - x*w);
            result[1][3] = T(0);

            result[2][0] = T(2) * (x*z - y*w);
            result[2][1] = T(2) * (y*z + x*w);
            result[2][2] = T(1) - T(2) * (x*x + y*y);
            result[2][3] = T(0);

            result[3][0] = T(0);
            result[3][1] = T(0);
            result[3][2] = T(0);
            result[3][3] = T(1);

            return result;
        }
    };

    // Quaternion specialization.
    typedef Tquaternion<float > quaternion;
    typedef Tquaternion<int> iquaternion;

    // Matrixes
    template <typename T, const int width, const int height> 
    class matXY
    {
    protected:
        // Arrays of a vector form a matrix.
        vecX<T, height> data[width];

    public:

        inline class vecX<T, height>& operator [] (int n)
        {
            return data[n];
        }
        inline const class vecX<T, height>& operator [] (int n) const
        {
            return data[n];
        }

        // Constructors
        inline matXY()
        { }
        inline matXY(const T clear)
        {
            for (int n = 0; n < width; n++)
            {
                for (int m = 0; m < height; m++)
                {
                    data[n][m] = clear;
                }
            }
        }
        inline matXY(const matXY& other)
        {
            for (int n = 0; n < width; n++)
            {
                data[n] = other.data[n];
            }
        }
        inline matXY& operator = (const class matXY<T, width, height>& other)
        {
            for (int n = 0; n < width; n++)
            {
                data[n] = other.data[n];
            }

            return *this;
        }

        // Data access functions
        inline operator T* () 
        {
            return &data[0][0];
        }
        inline operator const T* () const
        {
            return &data[0][0];
        }

        // Basic Mathematical operators

        // Matrix transpose
        inline matXY<T, height, width> transpose() const
        {
            matXY<T, height, width> result;

            for (int y = 0; y < width; y++)
            {
                for (int x = 0; x < height; x++)
                {
                    result[x][y] = data[y][x];
                }
            }

            return result;
        }
    
        // Identity matrix
        static inline class matXY<T, width, height> identity()
        {
            class matXY<T, width, height> result (0);

            for (int i = 0; i < width; i++)
            {
                result[i][i] = 1;
            }

            return result;
        }

        static inline int Width()
        {
            return width;
        }
        static inline int Height()
        {
            return height;
        }

        // Addition
        inline matXY& operator + (const class matXY<T, width, height>& other) const
        {
            class matXY<T, w, n> result;
            for (int n = 0; n < width; n++)
            {
                result[n] = data[n] + other[n];
            }

            return result;
        }
        inline matXY& operator += (const class matXY<T, width, height>& other)
        {
            for (int n = 0; n < width; n++)
            {
                data[n] += other[n];
            }

            return *this;
        }

        // Subtraction
        inline matXY& operator - (const class matXY<T, width, height>& other) const
        {
            class matXY<T, w, n> result;

            for (int n = 0; n < width; n++)
            {
                result[n] = data[n] - other[n];
            }

            return result;
        }
        inline matXY& operator -= (const class matXY<T, width, height>& other)
        {
            for (int n = 0; n < width; n++)
            {
                data[n] -= other[n];
            }

            return *this;
        }

        // Memberwise scalar multiplication
        inline matXY& operator * (const T other) const
        {
            class matXY<T, w, n> result;
            for (int n = 0; n < width; n++)
            {
                result[n] = data[n] * other;
            }

            return result;
        }
        inline matXY& operator += (const T other)
        {
            for (int n = 0; n < width; n++)
            {
                data[n] *= other;
            }

            return *this;
        }

        // Matrix multiplication (square matrixes only)
        inline class matXY<T, width, height> operator * (const class matXY<T, width, height>& other) const
        {
            class matXY<T, width, height> result (0);

            for (int n = 0; n < width; n++)
            {
                for (int m = 0; m < height; m++)
                {
                    for (int o = 0; o < width; o++)
                    {
                        result[n][m] += (data[o][m]*other[n][o]);
                    }
                }
            }

            return result;
        }
        inline class matXY<T, width, height>& operator *= (const class matXY<T, width, height>& other)
        {
            return (*this = *this * other);
        }
    };

    // Specialized 2x3 matrix.
    template <typename T>
    class Tmat2 : public matXY<T, 2, 2>
    {
    public: 
        inline Tmat2()
        { }
        inline Tmat2(const T clear) : matXY<T, 2, 2>(clear)
        { }
        inline Tmat2(const matXY<T, 2, 2>& other) : matXY<T, 2, 2>(other)
        { }
    };

    // Specialized 4x4 matrix.
    template <typename T>
    class Tmat4 : public matXY<T, 4, 4>
    {
    public: 
        inline Tmat4()
        { }
        inline Tmat4(const T clear) : matXY<T, 4, 4>(clear)
        { }
        inline Tmat4(const matXY<T, 4, 4>& other) : matXY<T,4, 4>(other)
        { }
    };

    // Definitions for float ing-point and integer matrixes
    typedef Tmat2<float > mat2;
    typedef Tmat2<int> imat2;
    
    typedef Tmat4<float > mat4;
    typedef Tmat4<int> imat4;

    // Construction of matrices for graphical display.

    // Generates a perspective frustrum projection matrix
    inline mat4 perspective(float  angleYDeg, float  aspect, float  nearP, float  farP)
    {
        mat4 result;

        result[0] = vec4(1.0f/(aspect*tan(ToRadians(0.5f*angleYDeg))), 0.0f, 0.0f, 0.0f);
        result[1] = vec4(0.0f, 1.0f/tan(ToRadians(0.5f*angleYDeg)), 0.0f, 0.0f);
        result[2] = vec4(0.0f, 0.0f, (nearP + farP)/(nearP - farP), -1.0f);
        result[3] = vec4(0.0f, 0.0f, (2.0f * nearP * farP) / (nearP - farP), 0.0f);

        return result;
    }

    // Generates an orthographic projection matrix
    inline mat4 orthographic(float  left, float  right, float  bottom, float  top, float  nearP, float  farP)
    {
        mat4 result;

        result[0] = vec4(2.0f / (right - left), 0.0f, 0.0f, 0.0f);
        result[1] = vec4(0.0f, 2.0f / (top - bottom), 0.0f, 0.0f);
        result[2] = vec4(0.0f, 0.0f, 2.0f / (nearP - farP), 0.0f);
        result[3] = vec4((left + right) / (left - right), (bottom + top) / (bottom - top), 
            (nearP + farP) / (farP - nearP), 1.0f);
        
        return result;
    }

    // Generates a look-at matrix
    inline mat4 lookat(const vec3 target, const vec3 camera, const vec3 up)
    {
        const vec3 forward = (target - camera).normalize();
        const vec3 upNorm = up.normalize();

        const vec3 side = forward.cross(upNorm);
        const vec3 upNew = side.cross(forward);

        mat4 result;

        result[0] = vec4(side[0], side[1], side[2], 0.0f);
        result[1] = vec4(upNew[0], upNew[1], upNew[2], 0.0f);
        result[2] = vec4(forward[0], forward[1], forward[2], 0.0f);
        result[3] = vec4(-camera[0], -camera[1], -camera[2], 1.0f);

        return result;
    }

    // Generates a translation matrix
    inline mat4 translate(const vec3 translation)
    {
        mat4 result;

        result[0] = vec4(1.0f, 0.0f, 0.0f, 0.0f);
        result[1] = vec4(0.0f, 1.0f, 0.0f, 0.0f);
        result[2] = vec4(0.0f, 0.0f, 1.0f, 0.0f);
        result[3] = vec4(translation[0], translation[1], translation[2], 1.0f);

        return result;
    }

    // Generates a scaling matrix from a vector
    inline mat4 scale(const vec3 scale)
    {
        mat4 result;

        result[0] = vec4(scale[0], 0.0f, 0.0f, 0.0f);
        result[1] = vec4(0.0f, scale[1], 0.0f, 0.0f);
        result[2] = vec4(0.0f, 0.0f, scale[2], 0.0f);
        result[3] = vec4(0.0f, 0.0f, 0.0f, 1.0f);

        return result;
    }

    // Generates a scaling matrix from a scalar
    inline mat4 scale(float  scale)
    {
        mat4 result;

        result[0] = vec4(scale, 0.0f, 0.0f, 0.0f);
        result[1] = vec4(0.0f, scale, 0.0f, 0.0f);
        result[2] = vec4(0.0f, 0.0f, scale, 0.0f);
        result[3] = vec4(0.0f, 0.0f, 0.0f, 1.0f);

        return result;
    }

    // Generates several forms of rotation matrixes, copied over from vmath.h, angle in degrees
    inline mat4 rotate(float  angle, float  x, float  y, float  z)
    {
        mat4 result;

        float  rads = ToRadians(angle);
        const float  c = cosf(rads);
        const float  s = sinf(rads);
        const float  omc = 1.0f - c;

        result[0] = vec4(( x * x * omc + c), (y * x * omc + z * s), (x * z * omc - y * s), 0.0f);
        result[1] = vec4((x * y * omc - z * s), (y * y * omc + c), (y * z * omc + x * s), 0.0f);
        result[2] = vec4((x * z * omc + y * s), (y * z * omc - x * s), (z* z * omc + c), 0.0f);
        result[3] = vec4(0.0f, 0.0f, 0.0f, 1.0f);

        return result;
    }
    inline mat4 rotate(float  angle, const vec3 axis)
    {
        return rotate(angle, axis[0], axis[1], axis[2]);
    }
    inline mat4 rotate(float  angleX, float  angleY, float  angleZ)
    {
        return rotate(angleZ, 0.0f, 0.0f, 1.0f) * rotate(angleY, 0.0f, 1.0f, 0.0f) * rotate(angleX, 1.0f, 0.0f, 0.0f);
    }
}
