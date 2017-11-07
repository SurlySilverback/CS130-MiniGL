/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>

using namespace std;

/**
 * Useful data types
 */
typedef mat<MGLfloat,4,4> mat4; //data structure storing a 4x4 matrix, see mat.h
typedef mat<MGLfloat,3,3> mat3; //data structure storing a 3x3 matrix, see mat.h
typedef vec<MGLfloat,4> vec4;   //data structure storing a 4 dimensional vector, see vec.h
typedef vec<MGLfloat,3> vec3;   //data structure storing a 3 dimensional vector, see vec.h
typedef vec<MGLfloat,2> vec2;   //data structure storing a 2 dimensional vector, see vec.h

/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}


/*************** Basic Structures ***************/

// Vertex
struct vertex
{
    /*Local Variables*/
    vec4 pos;
    vec3 color;

    /*Constructor with params*/
    vertex(vec4 pos, vec3 color) : pos(pos), color(color) {}
};


// Triangle
struct triangle
{
    /*Local Variables*/
    vertex a, b, c;
    
    /*Constructor with Params*/
    triangle(vertex v1, vertex v2, vertex v3) : a(v1), b(v2), c(v3) {}
};


/********** Variables **********/

// drawMode
MGLpoly_mode drawMode;

// listOfVertices
vector<vertex> listOfVertices;

// currColor
vec3 currColor;

// listOfTriangles
vector<triangle> listOfTriangles;


/********** Functions **********/

void mglColor( vec3 color )
{
    currColor = color;
};


/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{
};

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
    drawMode = mode;
};


/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
    int i = 0;
    int j = listOfVertices.size(); 

    switch( drawMode )
    {
        case MGL_TRIANGLES:

            // Take three vertices at a time from LoV and turn them into triangles in LoT
            for ( unsigned i = 0; i < listOfVertices.size(); i+=3 )
            {
                triangle newTriangle( listOfVertices.at( i ), listOfVertices.at( i+1 ), listOfVertices.at( i+2 ) );

                listOfTriangles.push_back( newTriangle );
            }

            listOfVertices.clear();

            break;

        case MGL_QUADS:

            // Take four vertices at a time from LoV and turn them into quads in LoT
            for ( unsigned i = 0; i < listOfVertices.size(); i+=4 )
            {
                // We need to handle cases
                



                /* Old Code (doesn't handle possible overlap)
                triangle newTriangle1( listOfVertices.at( i ), listOfVertices.at( i+1 ), listOfVertices.at( i+2 ) );
                triangle newTriangle2( listOfVertices.at( i+1 ), listOfVertices.at( i+2 ), listOfVertices.at( i+3 ) );
            
                listOfTriangles.push_back( newTriangle1 );
                listOfTriangles.push_back( newTriangle2 );
                */
            }

            listOfVertices.clear();

            break;
    }
};

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
    vec4 new_pos( x, y, 0, 1 );
    vertex v( new_pos,currColor );

    listOfVertices.push_back(v);
};

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
    vec4 new_pos( x, y, z, 1 );
    vertex v( new_pos, currColor );
    
    listOfVertices.push_back(v);
};

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
};

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
};

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
};

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
};

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
};

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
};

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
};

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
};

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
};

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
};

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
};

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
};
