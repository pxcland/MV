/*	-----------------------------------------------------------------------------
 *	MV - A Matrix Vector Mathematics Library
 *	Performs common and useful mathematical operations used in 3D graphics.
 *	C89 Compliant
 *	www.setsunasoft.com
 *
 *
 *	Copyright 2017 Patrick Cland
 *	
 *
 *	Permission is hereby granted, free of charge, to any person obtaining 
 *	a copy of this software and associated documentation files (the "Software"),
 *	to deal in the Software without restriction, including without limitation
 *	the rights to use, copy, modify, merge, publish, distribute, sublicense, 
 *	and/or sell copies of the Software, and to permit persons to whom the 
 *	Software is furnished to do so, subject to the following conditions:
 *
 *	The above copyright notice and this permission notice shall be included
 *	in all copies or substantial portions of the Software.

 *	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
 *	OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 *	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
 *	THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
 *	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
 *	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
 *	SOFTWARE.
 *	----------------------------------------------------------------------------- */

#ifndef MV_H
#define MV_H

#include <math.h>
#include <memory.h>
#include <stdio.h>

#ifdef __cplusplus
	extern "C" {
#endif



/* Vector Data Types */
typedef float MVvec2[2];
typedef float MVvec3[3];
typedef float MVvec4[4];

/* Matrix Data Types - Matrices are column major! */
typedef float MVmat2[4];
typedef float MVmat3[9];
typedef float MVmat4[16];

/* Useful constants */
#define MV_PI			(3.14159265359f)
#define MV_2PI			(6.28318530718f)
#define MV_PIDIV2		(1.57079632679f)
#define MV_PIDIV4		(0.78539816340f)
#define MV_PIDIV180		(0.01745329252f)
#define MV_PIDIV180I	(57.2957795131f)	/* Inverse of pi/180 */

/* Convert radians and degrees */
#define mvDegToRad(x)	((x) * MV_PIDIV180)
#define mvRadToDeg(x)	((x) * MV_PIDIV180I)


/* ------------------------------------------------------------------- */
/* ---------------------- Vector Operations -------------------------- */
/* ------------------------------------------------------------------- */

/* ATTENTION!! */
/* Functions which take a destination with a type of void* or float* */
/* don't protect against writing out of bounds. For example, copying */
/* an MVvec4 to an MVvec2 will overwrite 2 extra floats of whatever. */

/* set vector */
#define mvVecSet2(v, x, y)			v[0] = x; v[1] = y;
#define mvVecSet3(v, x, y, z)		v[0] = x; v[1] = y; v[2] = z;
#define mvVecSet4(v, x, y, z, w)	v[0] = x; v[1] = y; v[2] = z; v[3] = w;

/* copy vector */
void mvVecCopy2(void* dst, const void* src) { memcpy(dst, src, sizeof(MVvec2)); }
void mvVecCopy3(void* dst, const void* src) { memcpy(dst, src, sizeof(MVvec3)); }
void mvVecCopy4(void* dst, const void* src) { memcpy(dst, src, sizeof(MVvec4)); }

/* vector addition and subtraction: r = a +/- b */
void mvVecAdd2(float* r, const float* a, const float* b) { r[0] = a[0] + b[0]; r[1] = a[1] + b[1]; }
void mvVecAdd3(float* r, const float* a, const float* b) { r[0] = a[0] + b[0]; r[1] = a[1] + b[1]; r[2] = a[2] + b[2]; }
void mvVecAdd4(float* r, const float* a, const float* b) { r[0] = a[0] + b[0]; r[1] = a[1] + b[1]; r[2] = a[2] + b[2]; r[3] = a[3] + b[3]; }

void mvVecSub2(float* r, const float* a, const float* b) { r[0] = a[0] + b[0]; r[1] = a[1] + b[1]; }
void mvVecSub3(float* r, const float* a, const float* b) { r[0] = a[0] + b[0]; r[1] = a[1] + b[1]; r[2] = a[2] + b[2]; }
void mvVecSub4(float* r, const float* a, const float* b) { r[0] = a[0] + b[0]; r[1] = a[1] + b[1]; r[2] = a[2] + b[2]; r[3] = a[3] + b[3]; }

/* scale vector by factor */
void mvVecScale2(float* v, float factor) { v[0] *= factor; v[1] *= factor; }
void mvVecScale3(float* v, float factor) { v[0] *= factor; v[1] *= factor; v[2] *= factor;}
void mvVecScale4(float* v, float factor) { v[0] *= factor; v[1] *= factor; v[2] *= factor; v[3] *= factor; }

/* cross product */
void mvVecCrossProduct(float* r, const float* a, const float* b)
{
	r[0] = (a[1] * b[2]) - (b[1] * a[2]);
	r[1] = (-a[0] * b[2]) + (b[0] * a[2]);
	r[2] = (a[0] * b[1]) - (b[0] * a[1]);
}

/* scalar product */
float mvVecScalarProduct2(const float* a, const float* b) { return (a[0] * b[0]) + (a[1] * b[1]); }
float mvVecScalarProduct3(const float* a, const float* b) { return (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]); }
float mvVecScalarProduct4(const float* a, const float* b) { return (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]) + (a[3] * b[3]); }

/* vector length */
float mvVecLength2(const float* v) { return (float)sqrt((v[0] * v[0]) + (v[1] * v[1])); }
float mvVecLength3(const float* v) { return (float)sqrt((v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2])); }


/* angle between two vectors */
float mvVecAngleBetween2(const float* a, const float* b)
{
	return (float)acos( ((a[0] * b[0]) + (a[1] * b[1])) / ((float)sqrt((a[0] * a[0]) + (a[1] * a[1])) * sqrtf((b[0] * b[0]) + (b[1] * b[1]))));
}
float mvVecAngleBetween3(const float* a, const float* b)
{
	return (float)acos( ((a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2])) / ((float)sqrt((a[0] * a[0]) + (a[1] * a[1]) + (a[2] * a[2])) * sqrtf((b[0] * b[0]) + (b[1] * b[1]) + (b[2] * b[2]))));
}

/* normalize vector*/
void mvVecNormalize2(float* v)
{
	float length = (float)sqrt((v[0] * v[0]) + (v[1] * v[1]));
	v[0] /= length;
	v[1] /= length;
}
void mvVecNormalize3(float* v)
{
	float length = (float)sqrt((v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]));
	v[0] /= length;
	v[1] /= length;
	v[2] /= length;
}
void mvVecNormalize4(float* v)
{
	float length = (float)sqrt((v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]) + (v[3] * v[3]));
	v[0] /= length;
	v[1] /= length;
	v[2] /= length;
	v[3] /= length;
}

/* distance between vectors */
float mvVecDistanceBetween2(const float* a, const float* b)
{
	float tmp1 = a[0] - b[0];
	float tmp2 = a[1] - b[1];
	return (float)sqrt((tmp1 * tmp1) + (tmp2 * tmp2));
}
float mvVecDistanceBetween3(const float* a, const float* b)
{
	float tmp1 = a[0] - b[0];
	float tmp2 = a[1] - b[1];
	float tmp3 = a[2] - b[2];
	return (float)sqrt((tmp1 * tmp1) + (tmp2 * tmp2) + (tmp3 * tmp3));
}

/* get unit normal vector of a plane from 3 points that lie on the plane */
void mvVecUnitNormalFromPlane(void* r, const float* a, const float* b, const float* c)
{
	MVvec3 ab, ac;
	mvVecSub3(ab, b, a);
	mvVecSub3(ac, c, a);
	mvVecCrossProduct(r, ab, ac);
	mvVecNormalize3(r);
}



/* ------------------------------------------------------------------- */
/* ---------------------- Matrix Operations -------------------------- */
/* ------------------------------------------------------------------- */
/* copy matrix */
void mvMatCopy2(void* dst, const void* src) { memcpy(dst, src, sizeof(MVmat2)); }
void mvMatCopy3(void* dst, const void* src) { memcpy(dst, src, sizeof(MVmat3)); }
void mvMatCopy4(void* dst, const void* src) { memcpy(dst, src, sizeof(MVmat4)); }

/* load m with identity matrix */
void mvMatLoadIdentity2(MVmat2 m)
{
	MVmat2 tmp = {1.0f, 0.0f,
				  0.0f, 1.0f};
	memcpy(m, tmp, sizeof(MVmat2));
}
void mvMatLoadIdentity3(MVmat3 m)
{
	MVmat3 tmp = {1.0f, 0.0f, 0.0f,
				  0.0f, 1.0f, 0.0f,
				  0.0f, 0.0f, 1.0f};
	memcpy(m, tmp, sizeof(MVmat3));
}
void mvMatLoadIdentity4(MVmat4 m)
{
	MVmat4 tmp = {1.0f, 0.0f, 0.0f, 0.0f,
				  0.0f, 1.0f, 0.0f, 0.0f,
				  0.0f, 0.0f, 1.0f, 0.0f,
				  0.0f, 0.0f, 0.0f, 1.0f};
	memcpy(m, tmp, sizeof(MVmat4));
}

/* set a column of the matrix from vector src */
void mvMatSetCol2(MVmat2 m, const float* src, int col) { memcpy(m + (2*col), src, sizeof(MVvec2)); }
void mvMatSetCol3(MVmat3 m, const float* src, int col) { memcpy(m + (3*col), src, sizeof(MVvec3)); }
void mvMatSetCol4(MVmat4 m, const float* src, int col) { memcpy(m + (4*col), src, sizeof(MVvec4)); }

/* get column from matrix into vector */
void mvMatGetCol2(float* v, const float* src, int col) { memcpy(v, src + (2*col), sizeof(MVvec2)); }
void mvMatGetCol3(float* v, const float* src, int col) { memcpy(v, src + (3*col), sizeof(MVvec3)); }
void mvMatGetCol4(float* v, const float* src, int col) { memcpy(v, src + (4*col), sizeof(MVvec4)); }

/* get rotation matrix from 4x4 matrix */
void mvMatGetRotation(MVmat3 dst, const MVmat4 src)
{
	memcpy(dst, src, sizeof(float) * 3); /* X axis */
	memcpy(dst + 3, src + 4, sizeof(float) * 3); /* Y */
	memcpy(dst + 6, src + 8, sizeof(float) * 3); /* Z */
}

/* set rotation submatrix in a 4x4 matrix */
void mvMatSetRotation(MVmat4 dst, const MVmat3 src)
{
	memcpy(dst, src, sizeof(float) * 3); /* X axis */
	memcpy(dst + 4, src + 3, sizeof(float) * 3); /* Y */
	memcpy(dst + 8, src + 6, sizeof(float) * 3); /* Z */
}

/* Get translation vector from 4x4 matrix */
void mvMatGetTrans3(float* dst, const MVmat4 src) { memcpy(dst, src + 12, sizeof(MVvec3)); }
void mvMatGetTrans4(float* dst, const MVmat4 src) { memcpy(dst, src + 12, sizeof(MVvec4)); }

/* Set translation vector in 4x4 matrix */
void mvMatSetTrans3(MVmat4 dst, const float* src) { memcpy(dst + 12, src, sizeof(MVvec3)); }
void mvMatSetTrans4(MVmat4 dst, const float* src) { memcpy(dst + 12, src, sizeof(MVvec4)); }

/* Matrix multiplication */
void mvMatMultiply2(float* r, const float* a, const float* b)
{
	MVmat2 tmp;
	tmp[0] = a[0]*b[0] + a[2]*b[1];
	tmp[1] = a[1]*b[0] + a[3]*b[1];
	tmp[2] = a[0]*b[2] + a[2]*b[3];
	tmp[3] = a[1]*b[2] + a[3]*b[3];
	memcpy(r, tmp, sizeof(MVmat2));
}
void mvMatMultiply3(float* r, const float* a, const float* b)
{
	MVmat3 tmp;
	int i, j;
	for(i = 0; i < 3; i++)
		for(j = 0; j < 3; j++)
			tmp[(j*3)+i] = (a[0+i] * b[(j*3)]) + (a[3+i] * b[(j*3)+1]) + (a[6+i] * b[(j*3)+2]);
	memcpy(r, tmp, sizeof(MVmat3));
}
void mvMatMultiply4(float* r, const float* a, const float* b)
{
	MVmat4 tmp;
	int i, j;
	for(i = 0; i < 4; i++)
		for(j = 0; j < 4; j++)
			tmp[(j*4)+i] = (a[0+i] * b[(j*4)]) + (a[4+i] * b[(j*4)+1]) + (a[8+i] * b[(j*4)+2]) + (a[12+i] * b[(j*4)+3]);
	memcpy(r, tmp, sizeof(MVmat4));
}

/* scale matrix */
void mvMatScale2(float* r, float factor)
{
	int i;
	for(i = 0; i < 4; i++)
		r[i] *= factor;
}
void mvMatScale3(float* r, float factor)
{
	int i;
	for(i = 0; i < 9; i++)
		r[i] *= factor;
}
void mvMatScale4(float* r, float factor)
{
	int i;
	for(i = 0; i < 16; i++)
		r[i] *= factor;
}

/* Create scaling matrix */
void mvMatCreateScale3(MVmat3 m, float x, float y, float z)
{
	MVmat3 tmp = {0.0f};
	tmp[0] = x;
	tmp[4] = y;
	tmp[8] = z;
	memcpy(m, tmp, sizeof(MVmat3));
}
void mvMatCreateScale4(MVmat4 m, float x, float y, float z)
{
	MVmat4 tmp = {0.0f};
	tmp[0] = x;
	tmp[5] = y;
	tmp[10] = z;
	tmp[15] = 1.0f;
	memcpy(m, tmp, sizeof(MVmat4));
}

/* Create rotation matrix */
void mvMatCreateRotation3(MVmat3 m, float angle, float xAxis, float yAxis, float zAxis)
{
	float x, y, z;
	float c = (float)cos(angle);
	float s = (float)sinf(angle);
	float t = 1.0f - (float)cos(angle);
	/* Get normalized axis */
	MVvec3 axis = {xAxis, yAxis, zAxis};
	mvVecNormalize3(axis);
	x = axis[0]; y = axis[1]; z = axis[2];
	{
		MVmat3 tmp = {	  (x*x*t)+c,      (t*x*y)-(s*z),   (t*x*z)+(s*y) \
						, (t*x*y)+(s*z),  (t*y*y)+c,       (t*y*z)-(s*x) \
						, (t*x*z)-(s*y),  (t*y*z)+(s*x),   (t*z*z)+c };
		memcpy(m, tmp, sizeof(MVmat3));
	}
}
void mvMatCreateRotation4(MVmat4 m, float angle, float xAxis, float yAxis, float zAxis)
{
	float x, y, z;
	float c = (float)cos(angle);
	float s = (float)sin(angle);
	float t = 1.0f - (float)cos(angle);
	/* Get normalized axis */
	MVvec3 axis = {xAxis, yAxis, zAxis};
	mvVecNormalize3(axis);
	x = axis[0]; y = axis[1]; z = axis[2];
	{
		MVmat4 tmp = {	  (x*x*t)+c,      (t*x*y)-(s*z),   (t*x*z)+(s*y),   0.0f \
						, (t*x*y)+(s*z),  (t*y*y)+c,       (t*y*z)-(s*x),   0.0f \
						, (t*x*z)-(s*y),  (t*y*z)+(s*x),   (t*z*z)+c,       0.0f \
						, 0.0f,           0.0f,            0.0f,            1.0f};
		memcpy(m, tmp, sizeof(MVmat4));
	}
}

/* rotate matrix */
void mvRotate(MVmat4 m, float angle, float x, float y, float z)
{
	MVmat3 rot, tmp;
	mvMatGetRotation(tmp, m);
	mvMatCreateRotation3(rot, angle, x, y, z);
	mvMatMultiply3(tmp, tmp, rot);
	mvMatSetRotation(m, tmp);
}

/* translate matrix */
void mvTranslate(MVmat4 m, float x, float y, float z) { m[12] += x; m[13] += y; m[14] += z; }

/* scale matrix */
void mvScale(MVmat4 m, float x, float y, float z)
{
	MVmat4 tmp;
	mvMatCreateScale4(tmp, x, y, z);
	mvMatMultiply4(tmp, m, tmp);
	memcpy(m, tmp, sizeof(MVmat4));
}

/* multiply vector by matrix */
void mvMatVecMultiply2(float* r, const float* m, const float* v)
{
	MVvec2 tmp;
	memcpy(tmp, v, sizeof(MVvec2));
	tmp[0] = m[0]*v[0] + m[2]*v[1];
	tmp[1] = m[1]*v[0] + m[3]*v[1];
	memcpy(r, tmp, sizeof(MVvec2));
}
void mvMatVecMultiply3(float* r, const float* m, const float* v)
{
	MVvec3 tmp;
	memcpy(tmp, v, sizeof(MVvec3));
	tmp[0] = m[0]*v[0] + m[3]*v[1] + m[6]*v[2];
	tmp[1] = m[1]*v[0] + m[4]*v[1] + m[7]*v[2];
	tmp[2] = m[2]*v[0] + m[5]*v[1] + m[8]*v[2];
	memcpy(r, tmp, sizeof(MVvec3));
}
void mvMatVecMultiply4(float* r, const float* m, const float* v)
{
	MVvec4 tmp;
	memcpy(tmp, v, sizeof(MVvec4));
	tmp[0] = m[0]*v[0] + m[4]*v[1] + m[8]*v[2] + m[12]*v[3];
	tmp[1] = m[1]*v[0] + m[5]*v[1] + m[9]*v[2] + m[13]*v[3];
	tmp[2] = m[2]*v[0] + m[6]*v[1] + m[10]*v[2] + m[14]*v[3];
	tmp[3] = m[3]*v[0] + m[7]*v[1] + m[11]*v[2] + m[15]*v[3];
	memcpy(r, tmp, sizeof(MVvec4));
}

/* Generate orthographic projection matrix */
void mvMatCreateOrthographic(MVmat4 m, float left, float right, float bottom, float top, float near, float far)
{
	MVmat4 tmp = {0.0f};
	tmp[0] = 2.0f/(right-left);
	tmp[5] = 2.0f/(top-bottom);
	tmp[10] = -2.0f/(far-near);
	tmp[12] = -(right+left)/(right-left);
	tmp[13] = -(top+bottom)/(top-bottom);
	tmp[14] = -(far+near)/(far-near);
	tmp[15] = 1.0f;
	memcpy(m, tmp, sizeof(MVmat4));
}

/* Multiply specified matrix with an orthographic projection matrix */
void mvOrtho(MVmat4 m, float left, float right, float bottom, float top, float near, float far)
{
	MVmat4 tmp;
	mvMatCreateOrthographic(tmp, left, right, bottom, top, near, far);
	mvMatMultiply4(tmp, m, tmp);
	memcpy(m, tmp, sizeof(MVmat4));
}

/* create perspective projection matrix */
void mvMatCreatePerspective(MVmat4 m, float left, float right, float bottom, float top, float near, float far)
{
	MVmat4 tmp = {0.0f};
	tmp[0] = (2.0f*near)/(right-left);
	tmp[5] = (2.0f*near)/(top-bottom);
	tmp[8] = (right+left)/(right-left);
	tmp[9] = (top+bottom)/(top-bottom);
	tmp[10] = -(far+near)/(far-near);
	tmp[11] = -1.0f;
	tmp[14] = -(2.0f*far*near)/(far-near);
	memcpy(m, tmp, sizeof(MVmat4));
}

/* multiply current matrix by perspective projection matrix */
void mvFrustum(MVmat4 m, float left, float right, float bottom, float top, float near, float far)
{
	MVmat4 tmp;
	mvMatCreatePerspective(tmp, left, right, bottom, top, near, far);
	mvMatMultiply4(tmp, m, tmp);
	memcpy(m, tmp, sizeof(MVmat4));
}

/* multiply current matrix by perspective projection matrix - angle in radians */
void mvPerspective(MVmat4 m, float fovyRadians, float aspect, float near, float far)
{
	float y, x;
	y = near * (float)tan(fovyRadians / 2.0f);
	x = y * aspect;
	mvFrustum(m, -x, x, -y, y, near, far);
}

/* ------------------------------------------------------------------- */
/* ---------------------- Stack Operations --------------------------- */
/* ------------------------------------------------------------------- */
#define MV_STACK_SIZE 32
typedef struct _MVstack
{
	MVmat4 matrix[MV_STACK_SIZE];
	int top;
} MVstack;

void mvPushMatrix(MVstack* s, const MVmat4 m)
{
	s->top++;
	memcpy(&(s->matrix[MV_STACK_SIZE - s->top]), m, sizeof(MVmat4));
}

void mvPopMatrix(MVstack* s, MVmat4 m)
{
	memcpy(m, &(s->matrix[MV_STACK_SIZE - s->top]), sizeof(MVmat4));
	s->top--;
}

void mvPeekMatrix(MVstack* s, MVmat4 m)
{
	memcpy(m, &(s->matrix[MV_STACK_SIZE - s->top]), sizeof(MVmat4));
}

int mvGetCurrentStackDepth(MVstack* s)
{
	return s->top;
}


/* ------------------------------------------------------------------- */
/* ---------------------- Console Output ----------------------------- */
/* ------------------------------------------------------------------- */
/* print an n dimension vector */
void mvVecPrint(const float* v, int dim)
{
	int i;
	if(dim < 1 || dim > 4)
	{
		puts("mvVecPrint() : invalid dimension; must be 2, 3, or 4.");
		return;
	}
	printf("<");
	for(i = 0; i < dim; i++)
	{
		printf(" %.2f ", v[i]);
	}
	puts(">");
}

/* print an nxn matrix */
void mvMatPrint(const float* m, int dim)
{
	int i, j;
	if(dim < 2 || dim > 4)
	{
		puts("mvVecPrint() : invalid dimension; must be 3 or 4");
		return;
	}
	puts("--------------------------------");
	for(i = 0; i < dim; i++)
	{
		for(j = 0; j < dim; j++)
		{
			printf(" %.2f\t", m[(j*dim)+i]);
		}
		putchar('\n');
	}
	puts("--------------------------------");
}

#ifdef __cplusplus
}
#endif

#endif
