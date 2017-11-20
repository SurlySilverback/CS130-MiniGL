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


/*******************************************************************************
 BASIC STRUCTURES
*******************************************************************************/

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

/*******************************************************************************
 VARIABLES
*******************************************************************************/

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

/*******************************************************************************
 DATA STRUCTURES
 ******************************************************************************/

// Projection Stack
stack<mat4> currentProjMatrix;

// Model View Stack
stack<mat4> currentModViewMatrix;

/*******************************************************************************
 FUNCTIONS
*******************************************************************************/

stack<mat4> * getMatrixModeRef()
{
  stack<mat4> *value;

  if ( matrixMode == MGL_PROJECTION )
    value = &currentProjMatrix;

  else
    value = &currentModViewMatrix;

  return value;
}

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
    //cout << "mglReadPixels(): listOfTriangles.size() = " << listOfTriangles.size() << endl << endl;
    // For every triangle inside of listOfTriangles...
    for ( unsigned n = 0; n < listOfTriangles.size(); ++n )
    {
      // ...grab the positions of its vertices and store it in a local triangle obj...
      triangle curr_triangle = listOfTriangles.at(n);

/*
      cout << "Before the object-to-display conversion." << endl << endl

           << "Triangle ABC coordinates:" << endl
           << "Vertex A: (" << curr_triangle.a.pos[0] << "," << curr_triangle.a.pos[1] << "," << curr_triangle.a.pos[2] << "," << curr_triangle.a.pos[3] << ")" << endl
           << "Vertex B: (" << curr_triangle.b.pos[0] << "," << curr_triangle.b.pos[1] << "," << curr_triangle.b.pos[2] << "," << curr_triangle.b.pos[3] << ")" << endl
           << "Vertex C: (" << curr_triangle.c.pos[0] << "," << curr_triangle.c.pos[1] << "," << curr_triangle.c.pos[2] << "," << curr_triangle.c.pos[3] << ")" << endl << endl;
*/

      // ...translate the object coords of the temp vertices into display coords...
      curr_triangle.a.pos[0] = ( width / ( 2 * curr_triangle.a.pos[3] ) ) * ( curr_triangle.a.pos[0] + 1 );
      curr_triangle.a.pos[1] = ( height / ( 2 * curr_triangle.a.pos[3] ) ) * ( curr_triangle.a.pos[1] + 1 );
      curr_triangle.b.pos[0] = ( width / ( 2 * curr_triangle.b.pos[3] ) ) * ( curr_triangle.b.pos[0] + 1 );
      curr_triangle.b.pos[1] = ( height / ( 2 * curr_triangle.b.pos[3] ) ) * ( curr_triangle.b.pos[1] + 1 );
      curr_triangle.c.pos[0] = ( width / ( 2 * curr_triangle.b.pos[3] ) ) * ( curr_triangle.c.pos[0] + 1 );
      curr_triangle.c.pos[1] = ( height / ( 2 * curr_triangle.b.pos[3] ) ) * ( curr_triangle.c.pos[1] + 1 );

/*
      cout << "After the object-to-display conversion." << endl << endl
           << "Triangle ABC coordinates:" << endl
           << "Vertex A: (" << curr_triangle.a.pos[0] << "," << curr_triangle.a.pos[1] << "," << curr_triangle.a.pos[2] << "," << curr_triangle.a.pos[3] << ")" << endl
           << "Vertex B: (" << curr_triangle.b.pos[0] << "," << curr_triangle.b.pos[1] << "," << curr_triangle.b.pos[2] << "," << curr_triangle.b.pos[3] << ")" << endl
           << "Vertex C: (" << curr_triangle.c.pos[0] << "," << curr_triangle.c.pos[1] << "," << curr_triangle.c.pos[2] << "," << curr_triangle.c.pos[3] << ")" << endl << endl;
*/

      // ...and then transform those temp coords into the bounding box.
      // Note: We declare our bounding vars as integers for the control loop that
      // follows for barycentric calculations. We perform casts since the values
      // we are working with inside of the temp triangle are MGLfloats.
      MGLint xmin = (MGLint)floor( min( min( curr_triangle.a.pos[0], curr_triangle.b.pos[0] ), curr_triangle.c.pos[0] ) );
      MGLint xmax = (MGLint)ceil( max( max( curr_triangle.a.pos[0], curr_triangle.b.pos[0] ), curr_triangle.c.pos[0] ) );
      MGLint ymin = (MGLint)floor( min( min( curr_triangle.a.pos[1], curr_triangle.b.pos[1] ), curr_triangle.c.pos[1] ) );
      MGLint ymax = (MGLint)ceil( max( max( curr_triangle.a.pos[1], curr_triangle.b.pos[1] ), curr_triangle.c.pos[1] ) );

      cout << "After the bounding box calculation." << endl
           << "(xmin,xmax) = (" << xmin << "," << xmax << ")" << endl
           << "(ymin,ymax) = (" << ymin << "," << ymax << ")" << endl << endl;

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

      //cout << "After creating vectors AB, AC, and BC." << endl << endl;

      // We now have curr_triangle's area.
      curr_triangle_area = ( (cross(ab, ac)).magnitude() ) / 2;

      //cout << "After finding area of ABC." << endl << endl;

      // For the constraints of our bounding box...
      for ( int i = xmin; i < xmax; ++i )
      {
        for ( int j = ymin; j < ymax; ++j )
        {
          //cout << "Beginning for-loop for (i,j) = (" << i << ',' << j << ")" << endl;
          vec4 p_pos(i,j,0,1);
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

          //cout << "After calculating vectors AP and BP." << endl << endl;

          // Now, get the areas of subtriangles ABP, ACP, and BCP
          areaABP = ( (cross(ab, ap)).magnitude() ) / 2;
          areaACP = ( (cross(ac, ap)).magnitude() ) / 2;
          areaBCP = ( (cross(bc, bp)).magnitude() ) / 2;

          //cout << "After calculating subtriangle areas." << endl << endl;

          /* NOTE: ANOTHER WAY TO FIND AREA - POSSIBLY FASTER
             Consider using this to speed up raster for the stress test

          area(abc) = a_x * ( b_y - c_y ) + a_y * ( c_x - b_x ) + ( b_x * c_y - b_y * c_x );

          ********************************************/

          // Find barycentric coordinates.
          float alpha = areaBCP / curr_triangle_area;
          float beta = areaACP / curr_triangle_area;
          float gamma = areaABP / curr_triangle_area;

          //cout << "After calculating alpha, beta, gamma." << endl;

          if ( alpha + beta + gamma <= 1.015f )
          {
              //cout << "(width,height) are (" << width << ", " << height << ")" << endl;
              //cout << "Setting pixel (" << i << "," << j << ")" << endl;
              *(data + i + j * width) = Make_Pixel(255,255,255);
              //cout << "Pixel (" << i << "," << j << ") was set successfully." << endl << endl;
          }

          //cout << "Completed for-loop for (i,j) = (" << i << ',' << j << ")" << endl << endl;
        }
      }

      //cout << "After the readPixel double-for loop, at the end." << endl << endl;
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
                triangle newTriangle2( listOfVertices.at( i ), listOfVertices.at( i+2 ), listOfVertices.at( i+3 ) );

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
    // Create a container for the new 2D vertex.
    vec4 new_pos( x, y, 0, 1 );

    // If the model view matrix stack is not empty, multiply the vertex by the matrix
    // at the top of the stack.
    if ( !( currentModViewMatrix.empty() ) )
    {
      new_pos = currentModViewMatrix.top() * new_pos;
    }

    // If the projection matrix stack is not empty, multiply the vertex by the matrix
    // at the top of the stack.
    if ( !( currentProjMatrix.empty() ) )
    {
      new_pos = currentProjMatrix.top() * new_pos;
    }

    // Give the vertex colour.
    vertex v( new_pos,currColor );

    // Push it into the list of vertices.
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
    // Create a container for the new 2D vertex.
    vec4 new_pos( x, y, z, 1 );

    // If the model view matrix stack is not empty, multiply the vertex by the matrix
    // at the top of the stack.
    if ( !( currentModViewMatrix.empty() ) )
    {
      new_pos = currentModViewMatrix.top() * new_pos;
    }

    // If the projection matrix stack is not empty, multiply the vertex by the matrix
    // at the top of the stack.
    if ( !( currentProjMatrix.empty() ) )
    {
      new_pos = currentProjMatrix.top() * new_pos;
    }

    // Give the vertex colour.
    vertex v( new_pos, currColor );

    // Push it into the list of vertices.
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

  stack<mat4> *matRef = getMatrixModeRef();

  id.make_zero();

  id.values[0] = 1.0f;
  id.values[5] = 1.0f;
  id.values[10] = 1.0f;
  id.values[15] = 1.0f;

  //cout << "Inside mglLoadIdentity, printing matrix: " << endl << endl
  //     << id << endl << endl;

  if ( matRef->empty() )
  {
    matRef->push(id);
  }
  else
    matRef->top() = id;
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
  mat4 frustum;

  frustum.make_zero();

  frustum.values[0] = ( 2.0f * near ) / (right - left);
  frustum.values[5] = ( 2.0f * near ) / (top - bottom);
  frustum.values[8] = (right + left) / (right - left);
  frustum.values[9] = (top + bottom) / (top - bottom);
  frustum.values[10] = -(far + near) / (far - near);
  frustum.values[11] = -1.0f;
  frustum.values[14] = -( 2 * far * near ) / ( far - near );

  stack<mat4> *matRef = getMatrixModeRef();

  if ( matRef->empty() )
    mglLoadIdentity();

  matRef->top() = frustum * matRef->top();
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

  ortho.make_zero();

  ortho.values[0] = 2.0f / (right - left);
  ortho.values[5] = 2.0f / (top - bottom);
  ortho.values[10] = -(2.0f) / (far - near);
  ortho.values[12] = -(right + left) / (right - left);
  ortho.values[13] = -(top + bottom) / (top - bottom);
  ortho.values[14] = -(far + near) / (far - near);
  ortho.values[15] = 1.0f;

  stack<mat4> *matRef = getMatrixModeRef();

  if ( matRef->empty() )
    mglLoadIdentity();

  matRef->top() = ortho * matRef->top();
};

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
  currColor.values[0] = red; //red
  currColor.values[1] = green; //green
  currColor.values[2] = blue; //blue
};

void mglColor( vec3 color )
{
    currColor = color;
};
