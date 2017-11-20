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
#include <stack>

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

// matrixMode
MGLmatrix_mode matrixMode;

// listOfVertices
vector<vertex> listOfVertices;

// currColor
vec3 currColor;

// listOfTriangles
vector<triangle> listOfTriangles;


/********** Data Structures **********/

stack<mat4> currentMatrix;


/********** Functions **********/

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
    //cout << "listOfTriangles.size() = " << listOfTriangles.size() << endl << endl;
    // For every triangle inside of listOfTriangles...
    for ( unsigned n = 0; n < listOfTriangles.size(); ++n )
    {
      // ...grab the positions of its vertices and store it in a local triangle obj...
      triangle curr_triangle = listOfTriangles.at(n);

      // ...translate the object coords of the temp vertices into display coords...
      curr_triangle.a.pos[0] = ( width / 2 ) * ( curr_triangle.a.pos[0] + 1 );
      curr_triangle.a.pos[1] = ( height / 2 ) * ( curr_triangle.a.pos[1] + 1 );
      curr_triangle.b.pos[0] = ( width / 2 ) * ( curr_triangle.b.pos[0] + 1 );
      curr_triangle.b.pos[1] = ( height / 2 ) * ( curr_triangle.b.pos[1] + 1 );
      curr_triangle.c.pos[0] = ( width / 2 ) * ( curr_triangle.c.pos[0] + 1 );
      curr_triangle.c.pos[1] = ( height / 2 ) * ( curr_triangle.c.pos[1] + 1 );

      // ...and then transform those temp coords into the bounding box.
      // Note: We declare our bounding vars as integers for the control loop that
      // follows for barycentric calculations. We perform casts since the values
      // we are working with inside of the temp triangle are MGLfloats.
      MGLint xmin = (MGLint)floor( min( min( curr_triangle.a.pos[0], curr_triangle.b.pos[0] ), curr_triangle.c.pos[0] ) );
      MGLint xmax = (MGLint)ceil( max( max( curr_triangle.a.pos[0], curr_triangle.b.pos[0] ), curr_triangle.c.pos[0] ) );
      MGLint ymin = (MGLint)floor( min( min( curr_triangle.a.pos[1], curr_triangle.b.pos[1] ), curr_triangle.c.pos[1] ) );
      MGLint ymax = (MGLint)ceil( max( max( curr_triangle.a.pos[1], curr_triangle.b.pos[1] ), curr_triangle.c.pos[1] ) );

      // Now we pre-calculate the area of curr_triangle to help us with barycentric calculation.
      // Make two vectors out of (a,b) and (a,c), cross them, and divide by 2 to get curr_triangle's area.
      vec3 ab, ac, bc;
      float curr_triangle_area, areaABP, areaACP, areaBCP;
      vec3 p_color(255, 255, 255);

      // Create vector AB
      ab[0] = curr_triangle.b.pos[0] - curr_triangle.a.pos[0];
      ab[1] = curr_triangle.b.pos[1] - curr_triangle.a.pos[1];
      ab[2] = curr_triangle.b.pos[2] - curr_triangle.a.pos[2];

      // Create vector AC
      ac[0] = curr_triangle.c.pos[0] - curr_triangle.a.pos[0];
      ac[1] = curr_triangle.c.pos[1] - curr_triangle.a.pos[1];
      ac[2] = curr_triangle.c.pos[2] - curr_triangle.a.pos[2];

      // Create vector BC -- we'll use this later so make it here
      bc[0] = curr_triangle.c.pos[0] - curr_triangle.b.pos[0];
      bc[1] = curr_triangle.c.pos[1] - curr_triangle.b.pos[1];
      bc[2] = curr_triangle.c.pos[2] - curr_triangle.b.pos[2];

      // We now have curr_triangle's area.
      curr_triangle_area = ( (cross(ab, ac)).magnitude() ) / 2;

      //cout << "curr_triangle_area is " << curr_triangle_area << endl << endl;

      // For the constraints of our bounding box...
      for ( int i = xmin; i < xmax; ++i )
      {
        for ( int j = ymin; j < ymax; ++j )
        {
          vec4 p_pos(i,j,0,0);
          vertex p(p_pos, p_color);
          vec3 ap, bp;

          // Create vector AP
          ap[0] = p.pos[0] - curr_triangle.a.pos[0];
          ap[1] = p.pos[1] - curr_triangle.a.pos[1];
          ap[2] = p.pos[2] - curr_triangle.a.pos[2];

          // Create vector BP
          bp[0] = p.pos[0] - curr_triangle.b.pos[0];
          bp[1] = p.pos[1] - curr_triangle.b.pos[1];
          bp[2] = p.pos[2] - curr_triangle.b.pos[2];

          // Now, get the areas of subtriangles ABP, ACP, and BCP
          areaABP = ( (cross(ab, ap)).magnitude() ) / 2;
          areaACP = ( (cross(ac, ap)).magnitude() ) / 2;
          areaBCP = ( (cross(bc, bp)).magnitude() ) / 2;

          // Find barycentric coordinates.
          float alpha = areaBCP / curr_triangle_area;
          float beta = areaACP / curr_triangle_area;
          float gamma = areaABP / curr_triangle_area;

          // If the barycentric coordinates fall within parameters, draw the point.
          if ( alpha + beta + gamma <= 1.0 )
          {
              *(data + i + j * width) = Make_Pixel(255,255,255);
          }
        }
      }
    }
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
    //int i = 0;
    //int j = listOfVertices.size();

    switch( drawMode )
    {
        case MGL_TRIANGLES:

            // Take three vertices at a time from LoV and turn them into triangles in LoT
            for ( unsigned i = 0; (i + 2) < listOfVertices.size(); i+=3 )
            {
                triangle newTriangle( listOfVertices.at( i ), listOfVertices.at( i+1 ), listOfVertices.at( i+2 ) );

                listOfTriangles.push_back( newTriangle );
            }

            listOfVertices.clear();

            break;

        case MGL_QUADS:
            // Take four vertices at a time from LoV and turn them into quads in LoT
            for ( unsigned i = 0; (i + 3) < listOfVertices.size(); i+=4 )
            {
                // FIXME: This implementation does not handle the situation where the selected triangles
                // contain overlapping space within a quad. Compile and run this to get it working, but if the
                // final version has clipping issues, consider this logic as one of the solution spaces for
                // fixing the issue.
                triangle newTriangle1( listOfVertices.at( i ), listOfVertices.at( i+1 ), listOfVertices.at( i+2 ) );
                triangle newTriangle2( listOfVertices.at( i+1 ), listOfVertices.at( i+2 ), listOfVertices.at( i+3 ) );

                listOfTriangles.push_back( newTriangle1 );
                listOfTriangles.push_back( newTriangle2 );
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

    if ( !currentMatrix.empty() )
    {
      new_pos = currentMatrix.top() * new_pos;
    }

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

    if ( !currentMatrix.empty() )
    {
      new_pos = currentMatrix.top() * new_pos;
    }

    vertex v( new_pos, currColor );

    listOfVertices.push_back(v);
};

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
  matrixMode = mode;
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
  mat4 id;
  MGLfloat one = 1.0;

  id.make_zero();

  id.values[0] = one;
  id.values[5] = one;
  id.values[10] = one;
  id.values[15] = one;

  //cout << "id matrix is: " << endl
  //     << id << endl;

  currentMatrix.push(id);
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
  mat4 ortho;
  MGLfloat one = 1.0, two = 2.0;

  ortho.make_zero();

  ortho.values[0] = two / (right - left);
  ortho.values[5] = two / (top - bottom);
  ortho.values[10] = -(two) / (far - near);
  ortho.values[12] = -(right + left) / (right - left);
  ortho.values[13] = -(top + bottom) / (top - bottom);
  ortho.values[14] = -(far + near) / (far - near);
  ortho.values[15] = one;

  currentMatrix.top() = ortho * currentMatrix.top();
};

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{

  currColor.values[0] = red;
  currColor.values[1] = green;
  currColor.values[2] = blue;
};

void mglColor( vec3 color )
{
    currColor = color;
};
