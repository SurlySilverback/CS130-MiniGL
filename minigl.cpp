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

mat4 currentModelViewMatrix;

mat4 currentProjectionMatrix;

// Projection Stack
stack<mat4> projectionMatrixStack;

// Model View Stack
stack<mat4> modelViewMatrixStack;

// z_buffer to store the locations of drawn pixels in mglReadPixels.
// z_buffer is a 2D vector of MGLints, with default value of 2 indicating
// no pixel is drawn at that pixel location (since object space values are
// defined as -1 <= val <= 1).
vector< vector<MGLint> > z_buffer;

/*******************************************************************************
 FUNCTIONS
*******************************************************************************/

// area() takes in three vertex objects and returns the area of the triangle
MGLfloat area( vertex a, vertex b, vertex c)
{
    //area(abc) = a_x * ( b_y - c_y ) + a_y * ( c_x - b_x ) + ( b_x * c_y - b_y * c_x );
    return ( ( 0.5f * ( a.pos[0]*b.pos[1] - a.pos[1]*b.pos[0] + b.pos[0]*c.pos[1] - b.pos[1]*c.pos[0] + c.pos[0]*a.pos[1] - c.pos[1]*a.pos[0] ) ) );
}


// Returns a pointer to the current transformation matrix.
mat4 * getCurrentMatrix()
{
  mat4 *value;

  if ( matrixMode == MGL_PROJECTION )
    value = &currentProjectionMatrix;

  else
    value = &currentModelViewMatrix;

  return value;
}

// Returns a pointer to the current matrix view mode. All push/pop operations
// will be performed on this stack until the mode is changed using mglBegin()
stack<mat4> * getMatrixModeRef()
{
  stack<mat4> *value;

  if ( matrixMode == MGL_PROJECTION )
    value = &projectionMatrixStack;

  else
    value = &modelViewMatrixStack;

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
  cout << "mglReadPixels() begins here..." << endl;

  // Initialize the z-buffer to a default max for this draw pass.
  z_buffer.resize(width);

  for ( unsigned i = 0; i < width; ++i )
  {
    z_buffer.at(i).resize(height);

    for ( unsigned j = 0; j < height; ++j )
      z_buffer.at(i).at(j) = 2.0f;
  }
/*
    cout << "listOfTriangles.size() == " << listOfTriangles.size() << endl
         << "Value of *(data) is " << *(data) << endl << endl;
*/
    // For every triangle inside of listOfTriangles...
  for ( unsigned n = 0; n < listOfTriangles.size(); ++n )
  {
    // ...grab the positions of its vertices and store it in a local triangle obj...
    triangle curr_triangle = listOfTriangles.at(n);


    // Divide by the w value for homogeneous coordinates ( (x/w), (y/w), (z/w), 1 )
    // w value is stored at curr_triangle.a.pos[3]
    curr_triangle.a.pos = curr_triangle.a.pos / curr_triangle.a.pos[3];
    curr_triangle.b.pos = curr_triangle.b.pos / curr_triangle.b.pos[3];
    curr_triangle.c.pos = curr_triangle.c.pos / curr_triangle.c.pos[3];


    // ...translate the object coords of the temp vertices into display coords...
    // x = ( width / 2 ) ( x + 1 )
    // y = ( height / 2 ) ( y + 1 )
    curr_triangle.a.pos[0] = ( width * (0.5f) ) * ( curr_triangle.a.pos[0] + 1 );
    curr_triangle.a.pos[1] = ( height * (0.5f) ) * ( curr_triangle.a.pos[1] + 1 );
    curr_triangle.b.pos[0] = ( width * (0.5f) ) * ( curr_triangle.b.pos[0] + 1 );
    curr_triangle.b.pos[1] = ( height * (0.5f) ) * ( curr_triangle.b.pos[1] + 1 );
    curr_triangle.c.pos[0] = ( width * (0.5f) ) * ( curr_triangle.c.pos[0] + 1 );
    curr_triangle.c.pos[1] = ( height * (0.5f) ) * ( curr_triangle.c.pos[1] + 1 );

      cout << "After obj-to-display conversion:" << endl
           << "Triangle # " << n+1 << endl
           << "============" << endl
           << "A: (" << curr_triangle.a.pos[0] << ", " << curr_triangle.a.pos[1] << ", " << curr_triangle.a.pos[2] << ")" << endl
           << "B: (" << curr_triangle.b.pos[0] << ", " << curr_triangle.b.pos[1] << ", " << curr_triangle.b.pos[2] << ")" << endl
           << "C: (" << curr_triangle.c.pos[0] << ", " << curr_triangle.c.pos[1] << ", " << curr_triangle.c.pos[2] << ")" << endl << endl;

    // ...and then transform those temp coords into the bounding box.
    // Note: We declare our bounding vars as integers for the control loop that
    // follows for barycentric calculations. We perform casts since the values
    // we are working with inside of the temp triangle are MGLfloats.
    MGLint xmin = (MGLint)floor( min( min( curr_triangle.a.pos[0], curr_triangle.b.pos[0] ), curr_triangle.c.pos[0] ) );
    MGLint xmax = (MGLint)ceil( max( max( curr_triangle.a.pos[0], curr_triangle.b.pos[0] ), curr_triangle.c.pos[0] ) );
    MGLint ymin = (MGLint)floor( min( min( curr_triangle.a.pos[1], curr_triangle.b.pos[1] ), curr_triangle.c.pos[1] ) );
    MGLint ymax = (MGLint)ceil( max( max( curr_triangle.a.pos[1], curr_triangle.b.pos[1] ), curr_triangle.c.pos[1] ) );


    // Now we pre-calculate the area of curr_triangle to help us with barycentric calculation.
    float curr_triangle_area, areaABP, areaAPC, areaPBC;
    curr_triangle_area = area( curr_triangle.a, curr_triangle.b, curr_triangle.c );


    // FIXME: We will need to fix this hack for colour interpolation
    vec3 p_color = curr_triangle.a.color;


    // For the constraints of our bounding box...
    for ( int i = xmin; i < xmax; ++i )
    {
      for ( int j = ymin; j < ymax; ++j )
      {
        // FIXME: Here we must pass the interpolated value of z into the z-value for P
        vec4 p_pos(i,j,0,1);
        vertex p(p_pos, p_color);


        // Find the area of three component triangles inside ABC using point P at (i,j)
        areaABP = area( curr_triangle.a, curr_triangle.b, p );
        areaAPC = area( curr_triangle.a, p, curr_triangle.c );
        areaPBC = area( p, curr_triangle.b, curr_triangle.c );

        float alpha = areaPBC / curr_triangle_area;
        float beta = areaAPC / curr_triangle_area;
        float gamma = areaABP / curr_triangle_area;

// BEGIN: COLOUR INTERPOLATION BLOCK
        //float alpha_prime = alpha * curr_triangle.a.pos[2] *

// END: COLOUR INTERPOLATION BLOCK

        // Calculate current Z value for pixel (i,j)
        MGLfloat currZ = ( alpha*curr_triangle.a.pos[2] ) + ( beta*curr_triangle.b.pos[2] ) + ( gamma*curr_triangle.c.pos[2] );

        // Calculate the linear interpolation for colour
        vec3 interpColor;
        interpColor[0] = ( alpha*curr_triangle.a.color[0] ) + ( beta*curr_triangle.b.color[0] ) + ( gamma*curr_triangle.c.color[0] );
        interpColor[1] = ( alpha*curr_triangle.a.color[1] ) + ( beta*curr_triangle.b.color[1] ) + ( gamma*curr_triangle.c.color[1] );
        interpColor[2] = ( alpha*curr_triangle.a.color[2] ) + ( beta*curr_triangle.b.color[2] ) + ( gamma*curr_triangle.c.color[2] );

        // If all barycentric params are >=, then point p at (i,j) lies inside
        // the main triangle.
        if ( alpha >= 0 && beta >= 0 && gamma >= 0 )
        {
          // FIXME: Now check if pixel (i,j) has been drawn to. If it has, only draw
          // the current incoming pixel if its z-value is less than the one
          // already drawn at (i,j).
          if ( currZ < z_buffer.at(i).at(j) )
          {
            // Set the z_buffer equal to the new, closer pixel at (i,j)
            z_buffer.at(i).at(j) = currZ;

            // FIXME: Ensure that correct colour interpolation occurs here by passing params to Make_Pixel
            //*(data + i + j * width) = Make_Pixel( p_color[0] * 255, p_color[1] * 255, p_color[2] * 255 ); // OLD CODE
            *(data + i + j * width) = Make_Pixel( interpColor[0] * 255, interpColor[1] * 255, interpColor[2] * 255 );
          }
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
            cout << "End of mglEnd() in TRIANGLES" << endl << endl;

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
            cout << "End of mglEnd() in QUADS" << endl << endl;

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
  vec4 new_pos( x, y, 0.0f, 1.0f );


  cout << "In mglVertex2() -- currentModelViewMatrix is: "
       << currentModelViewMatrix << endl << endl;

  cout << "vertex coordinates before modView" << endl
       << "=================================" << endl
       << "x: " << new_pos[0] << endl
       << "y: " << new_pos[1] << endl
       << "z: " << new_pos[2] << endl
       << "w: " << new_pos[3] << endl << endl;

  // Multiply the new position by the currentModelViewMatrix
  new_pos = currentModelViewMatrix * new_pos;

  cout << "vertex coordinates after modView" << endl
       << "================================" << endl
       << "x: " << new_pos[0] << endl
       << "y: " << new_pos[1] << endl
       << "z: " << new_pos[2] << endl
       << "w: " << new_pos[3] << endl << endl;

  // Multiply the new position by the currentProjectionMatrix
  new_pos = currentProjectionMatrix * new_pos;

  cout << "vertex coordinates after projMatrix" << endl
       << "===================================" << endl
       << "x: " << new_pos[0] << endl
       << "y: " << new_pos[1] << endl
       << "z: " << new_pos[2] << endl
       << "w: " << new_pos[3] << endl << endl;

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
  // Create a position container for the new 3D vertex.
  vec4 new_pos( x, y, z, 1.0f );

  cout << "In mglVertex3() -- currentModelViewMatrix is: "
       << currentModelViewMatrix << endl << endl;

  cout << "vertex coordinates before modView" << endl
       << "=================================" << endl
       << "x: " << new_pos[0] << endl
       << "y: " << new_pos[1] << endl
       << "z: " << new_pos[2] << endl
       << "w: " << new_pos[3] << endl << endl;

  // Multiply the new position by the currentModelViewMatrix
  new_pos = currentModelViewMatrix * new_pos;

  cout << "vertex coordinates after modView" << endl
       << "================================" << endl
       << "x: " << new_pos[0] << endl
       << "y: " << new_pos[1] << endl
       << "z: " << new_pos[2] << endl
       << "w: " << new_pos[3] << endl << endl;

  // Multiply the new position by the currentProjectionMatrix
  new_pos = currentProjectionMatrix * new_pos;

  cout << "vertex coordinates after projMatrix" << endl
       << "===================================" << endl
       << "x: " << new_pos[0] << endl
       << "y: " << new_pos[1] << endl
       << "z: " << new_pos[2] << endl
       << "w: " << new_pos[3] << endl << endl;

  // Create the vertex using the adjusted position and the current colour.
  vertex v( new_pos, currColor );

  // Push the vertex into the list of vertices.
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
  stack<mat4> *matStackRef = getMatrixModeRef();
  mat4 *currMatRef = getCurrentMatrix();

  // FIXME: Make sure that popping this doesn't lead to a memory leak. Consider
  // invoking delete on pop operation in mglPopMatrix().
  matStackRef->push( *(currMatRef) );
};

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
  // FIXME: in mglPushMatrix, a dereferenced pointer is used to push matrices into
  // the current stack. Check here to make sure that the matrix being popped doesn't
  // also need to be deleted in order to avoid a possible memory leak.
  stack<mat4> *matStackRef = getMatrixModeRef();
  mat4 *currMatRef = getCurrentMatrix();

  if ( !( matStackRef->empty() ) )
  {
    // FIXME: This could be wrong. Originally, we just popped the top matrix of matStackRef without assignment.
    *(currMatRef) = matStackRef->top();
    matStackRef->pop();
  }
};

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
  // Create an identity matrix
  mat4 id = {{ 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f }};

  // Get the current matrix
  mat4 *currMatRef = getCurrentMatrix();

  *(currMatRef) = id;
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
  //currentMatrix = *(matrix);
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
  mat4 multMatrix;
  multMatrix.make_zero();

  mat4 *currMatRef = getCurrentMatrix();

  for (unsigned i = 0; i < 16; ++i)
  {
      multMatrix.values[i] = matrix[i];
  }

  *(currMatRef) = (*(currMatRef)) * multMatrix;
};

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
  mat4 translate = {{ 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, x, y, z, 1.0f }};

  mat4 *currMatRef = getCurrentMatrix();

  *(currMatRef) = (*(currMatRef)) * translate;
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
  mat4 rotate;
  MGLfloat x_ = x, y_ = y, z_ = z;

  // Check for normalization of x,y,z
  vec3 normTest(x_,y_,z_);

  // If the vector (x_, y_, z_) is not normalized, normalize it
  if ( normTest.magnitude() != 1 )
  {
    normTest.normalized();

    x_ = x_ / normTest.magnitude();
    y_ = y_ / normTest.magnitude();
    z_ = z_ / normTest.magnitude();
  }

  // Convert from degrees to radians before calculating sin, cos
  angle = ( angle * 3.14159265 ) / 180;

  float s = sin(angle);
  float c = cos(angle);

  rotate.make_zero();

  rotate.values[0] = pow(x_,2) * ( 1-c ) + c;
  rotate.values[1] = y_ * x_ * ( 1-c ) + z_ * s;
  rotate.values[2] = x_ * z_ * ( 1-c ) - y_ * s;
  rotate.values[4] = x_ * y_ * ( 1-c ) - z_ * s;
  rotate.values[5] = pow(y_,2) * ( 1-c ) + c;
  rotate.values[6] = y_ * z_ * ( 1-c ) + x_ * s;
  rotate.values[8] = x_ * z_ * ( 1-c ) + y_ * s;
  rotate.values[9] = y_ * z_ * ( 1-c ) - x_ * s;
  rotate.values[10] = pow(z_,2) * ( 1-c ) + c;
  rotate.values[15] = 1.0f;

  mat4 *currMatRef = getCurrentMatrix();

  *(currMatRef) = (*(currMatRef)) * rotate;
};

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
  mat4 scale = {{ x, 0.0f, 0.0f, 0.0f, 0.0f, y, 0.0f, 0.0f, 0.0f, 0.0f, z, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f }};

  mat4 *currMatRef = getCurrentMatrix();

  *(currMatRef) = (*(currMatRef)) * scale;
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

  // Create the perspective projection (frustum) matrix
  frustum.values[0] = ( 2.0f * near ) / ( right - left );
  frustum.values[5] = ( 2.0f * near ) / ( top - bottom );
  frustum.values[8] = ( right + left ) / ( right - left );
  frustum.values[9] = ( top + bottom ) / ( top - bottom );
  frustum.values[10] = -( far + near ) / ( far - near );
  frustum.values[11] = -1.0f;
  frustum.values[14] = -( 2.0f * far * near ) / ( far - near );

  mat4 *currMatRef = getCurrentMatrix();

  *(currMatRef) = (*(currMatRef)) * frustum;
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

  // Create the orthographic matrix
  ortho.values[0] = 2.0f / (right - left);
  ortho.values[5] = 2.0f / (top - bottom);
  ortho.values[10] = -(2.0f) / (far - near);
  ortho.values[12] = -(right + left) / (right - left);
  ortho.values[13] = -(top + bottom) / (top - bottom);
  ortho.values[14] = -(far + near) / (far - near);
  ortho.values[15] = 1.0f;

  mat4 *currMatRef = getCurrentMatrix();

  *(currMatRef) = (*(currMatRef)) * ortho;
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
