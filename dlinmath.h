#ifndef LINMATH_H
#define LINMATH_H

#include <dmath.h>

#define LINMATH_H_DEFINE_DVEC(n) \
typedef double dvec##n[n]; \
static inline void dvec##n##_add(dvec##n r, dvec##n const a, dvec##n const b) \
{ \
	int i; \
	for(i=0; i<n; ++i) \
		r[i] = a[i] + b[i]; \
} \
static inline void dvec##n##_sub(dvec##n r, dvec##n const a, dvec##n const b) \
{ \
	int i; \
	for(i=0; i<n; ++i) \
		r[i] = a[i] - b[i]; \
} \
static inline void dvec##n##_scale(dvec##n r, dvec##n const v, double const s) \
{ \
	int i; \
	for(i=0; i<n; ++i) \
		r[i] = v[i] * s; \
} \
static inline double dvec##n##_mul_inner(dvec##n const a, dvec##n const b) \
{ \
	double p = 0.; \
	int i; \
	for(i=0; i<n; ++i) \
		p += b[i]*a[i]; \
	return p; \
} \
static inline double dvec##n##_len(dvec##n const v) \
{ \
	return sqrtf(dvec##n##_mul_inner(v,v)); \
} \
static inline void dvec##n##_norm(dvec##n r, dvec##n const v) \
{ \
	double k = 1.0 / dvec##n##_len(v); \
	dvec##n##_scale(r, v, k); \
} \
static inline void dvec##n##_min(dvec##n r, dvec##n a, dvec##n b) \
{ \
	int i; \
	for(i=0; i<n; ++i) \
		r[i] = a[i]<b[i] ? a[i] : b[i]; \
} \
static inline void dvec##n##_max(dvec##n r, dvec##n a, dvec##n b) \
{ \
	int i; \
	for(i=0; i<n; ++i) \
		r[i] = a[i]>b[i] ? a[i] : b[i]; \
}

LINMATH_H_DEFINE_VEC(2)
LINMATH_H_DEFINE_VEC(3)
LINMATH_H_DEFINE_VEC(4)

static inline void dvec3_mul_cross(dvec3 r, dvec3 const a, dvec3 const b)
{
	r[0] = a[1]*b[2] - a[2]*b[1];
	r[1] = a[2]*b[0] - a[0]*b[2];
	r[2] = a[0]*b[1] - a[1]*b[0];
}

static inline void dvec3_reflect(dvec3 r, dvec3 const v, dvec3 const n)
{
	double p  = 2.f*dvec3_mul_inner(v, n);
	int i;
	for(i=0;i<3;++i)
		r[i] = v[i] - p*n[i];
}

static inline void dvec4_mul_cross(dvec4 r, dvec4 a, dvec4 b)
{
	r[0] = a[1]*b[2] - a[2]*b[1];
	r[1] = a[2]*b[0] - a[0]*b[2];
	r[2] = a[0]*b[1] - a[1]*b[0];
	r[3] = 1.f;
}

static inline void dvec4_reflect(dvec4 r, dvec4 v, dvec4 n)
{
	double p  = 2.f*dvec4_mul_inner(v, n);
	int i;
	for(i=0;i<4;++i)
		r[i] = v[i] - p*n[i];
}

typedef dvec4 dmat4x4[4];
static inline void dmat4x4_identity(dmat4x4 M)
{
	int i, j;
	for(i=0; i<4; ++i)
		for(j=0; j<4; ++j)
			M[i][j] = i==j ? 1.f : 0.f;
}
static inline void dmat4x4_dup(dmat4x4 M, dmat4x4 N)
{
	int i, j;
	for(i=0; i<4; ++i)
		for(j=0; j<4; ++j)
			M[i][j] = N[i][j];
}
static inline void dmat4x4_row(dvec4 r, dmat4x4 M, int i)
{
	int k;
	for(k=0; k<4; ++k)
		r[k] = M[k][i];
}
static inline void dmat4x4_col(dvec4 r, dmat4x4 M, int i)
{
	int k;
	for(k=0; k<4; ++k)
		r[k] = M[i][k];
}
static inline void dmat4x4_transpose(dmat4x4 M, dmat4x4 N)
{
	int i, j;
	for(j=0; j<4; ++j)
		for(i=0; i<4; ++i)
			M[i][j] = N[j][i];
}
static inline void dmat4x4_add(dmat4x4 M, dmat4x4 a, dmat4x4 b)
{
	int i;
	for(i=0; i<4; ++i)
		dvec4_add(M[i], a[i], b[i]);
}
static inline void dmat4x4_sub(dmat4x4 M, dmat4x4 a, dmat4x4 b)
{
	int i;
	for(i=0; i<4; ++i)
		dvec4_sub(M[i], a[i], b[i]);
}
static inline void dmat4x4_scale(dmat4x4 M, dmat4x4 a, double k)
{
	int i;
	for(i=0; i<4; ++i)
		dvec4_scale(M[i], a[i], k);
}
static inline void dmat4x4_scale_aniso(dmat4x4 M, dmat4x4 a, double x, double y, double z)
{
	int i;
	dvec4_scale(M[0], a[0], x);
	dvec4_scale(M[1], a[1], y);
	dvec4_scale(M[2], a[2], z);
	for(i = 0; i < 4; ++i) {
		M[3][i] = a[3][i];
	}
}
static inline void dmat4x4_mul(dmat4x4 M, dmat4x4 a, dmat4x4 b)
{
	dmat4x4 temp;
	int k, r, c;
	for(c=0; c<4; ++c) for(r=0; r<4; ++r) {
		temp[c][r] = 0.f;
		for(k=0; k<4; ++k)
			temp[c][r] += a[k][r] * b[c][k];
	}
	dmat4x4_dup(M, temp);
}
static inline void dmat4x4_mul_dvec4(dvec4 r, dmat4x4 M, dvec4 v)
{
	int i, j;
	for(j=0; j<4; ++j) {
		r[j] = 0.f;
		for(i=0; i<4; ++i)
			r[j] += M[i][j] * v[i];
	}
}
static inline void dmat4x4_translate(dmat4x4 T, double x, double y, double z)
{
	dmat4x4_identity(T);
	T[3][0] = x;
	T[3][1] = y;
	T[3][2] = z;
}
static inline void dmat4x4_translate_in_place(dmat4x4 M, double x, double y, double z)
{
	dvec4 t = {x, y, z, 0};
	dvec4 r;
	int i;
	for (i = 0; i < 4; ++i) {
		dmat4x4_row(r, M, i);
		M[3][i] += dvec4_mul_inner(r, t);
	}
}
static inline void dmat4x4_from_dvec3_mul_outer(dmat4x4 M, dvec3 a, dvec3 b)
{
	int i, j;
	for(i=0; i<4; ++i) for(j=0; j<4; ++j)
		M[i][j] = i<3 && j<3 ? a[i] * b[j] : 0.f;
}
static inline void dmat4x4_rotate(dmat4x4 R, dmat4x4 M, double x, double y, double z, double angle)
{
	double s = sinf(angle);
	double c = cosf(angle);
	dvec3 u = {x, y, z};

	if(dvec3_len(u) > 1e-4) {
		dvec3_norm(u, u);
		dmat4x4 T;
		dmat4x4_from_dvec3_mul_outer(T, u, u);

		dmat4x4 S = {
			{    0,  u[2], -u[1], 0},
			{-u[2],     0,  u[0], 0},
			{ u[1], -u[0],     0, 0},
			{    0,     0,     0, 0}
		};
		dmat4x4_scale(S, S, s);

		dmat4x4 C;
		dmat4x4_identity(C);
		dmat4x4_sub(C, C, T);

		dmat4x4_scale(C, C, c);

		dmat4x4_add(T, T, C);
		dmat4x4_add(T, T, S);

		T[3][3] = 1.;		
		dmat4x4_mul(R, M, T);
	} else {
		dmat4x4_dup(R, M);
	}
}
static inline void dmat4x4_rotate_X(dmat4x4 Q, dmat4x4 M, double angle)
{
	double s = sinf(angle);
	double c = cosf(angle);
	dmat4x4 R = {
		{1.f, 0.f, 0.f, 0.f},
		{0.f,   c,   s, 0.f},
		{0.f,  -s,   c, 0.f},
		{0.f, 0.f, 0.f, 1.f}
	};
	dmat4x4_mul(Q, M, R);
}
static inline void dmat4x4_rotate_Y(dmat4x4 Q, dmat4x4 M, double angle)
{
	double s = sinf(angle);
	double c = cosf(angle);
	dmat4x4 R = {
		{   c, 0.f,   s, 0.f},
		{ 0.f, 1.f, 0.f, 0.f},
		{  -s, 0.f,   c, 0.f},
		{ 0.f, 0.f, 0.f, 1.f}
	};
	dmat4x4_mul(Q, M, R);
}
static inline void dmat4x4_rotate_Z(dmat4x4 Q, dmat4x4 M, double angle)
{
	double s = sinf(angle);
	double c = cosf(angle);
	dmat4x4 R = {
		{   c,   s, 0.f, 0.f},
		{  -s,   c, 0.f, 0.f},
		{ 0.f, 0.f, 1.f, 0.f},
		{ 0.f, 0.f, 0.f, 1.f}
	};
	dmat4x4_mul(Q, M, R);
}
static inline void dmat4x4_invert(dmat4x4 T, dmat4x4 M)
{
	double s[6];
	double c[6];
	s[0] = M[0][0]*M[1][1] - M[1][0]*M[0][1];
	s[1] = M[0][0]*M[1][2] - M[1][0]*M[0][2];
	s[2] = M[0][0]*M[1][3] - M[1][0]*M[0][3];
	s[3] = M[0][1]*M[1][2] - M[1][1]*M[0][2];
	s[4] = M[0][1]*M[1][3] - M[1][1]*M[0][3];
	s[5] = M[0][2]*M[1][3] - M[1][2]*M[0][3];

	c[0] = M[2][0]*M[3][1] - M[3][0]*M[2][1];
	c[1] = M[2][0]*M[3][2] - M[3][0]*M[2][2];
	c[2] = M[2][0]*M[3][3] - M[3][0]*M[2][3];
	c[3] = M[2][1]*M[3][2] - M[3][1]*M[2][2];
	c[4] = M[2][1]*M[3][3] - M[3][1]*M[2][3];
	c[5] = M[2][2]*M[3][3] - M[3][2]*M[2][3];
	
	/* Assumes it is invertible */
	double idet = 1.0f/( s[0]*c[5]-s[1]*c[4]+s[2]*c[3]+s[3]*c[2]-s[4]*c[1]+s[5]*c[0] );
	
	T[0][0] = ( M[1][1] * c[5] - M[1][2] * c[4] + M[1][3] * c[3]) * idet;
	T[0][1] = (-M[0][1] * c[5] + M[0][2] * c[4] - M[0][3] * c[3]) * idet;
	T[0][2] = ( M[3][1] * s[5] - M[3][2] * s[4] + M[3][3] * s[3]) * idet;
	T[0][3] = (-M[2][1] * s[5] + M[2][2] * s[4] - M[2][3] * s[3]) * idet;

	T[1][0] = (-M[1][0] * c[5] + M[1][2] * c[2] - M[1][3] * c[1]) * idet;
	T[1][1] = ( M[0][0] * c[5] - M[0][2] * c[2] + M[0][3] * c[1]) * idet;
	T[1][2] = (-M[3][0] * s[5] + M[3][2] * s[2] - M[3][3] * s[1]) * idet;
	T[1][3] = ( M[2][0] * s[5] - M[2][2] * s[2] + M[2][3] * s[1]) * idet;

	T[2][0] = ( M[1][0] * c[4] - M[1][1] * c[2] + M[1][3] * c[0]) * idet;
	T[2][1] = (-M[0][0] * c[4] + M[0][1] * c[2] - M[0][3] * c[0]) * idet;
	T[2][2] = ( M[3][0] * s[4] - M[3][1] * s[2] + M[3][3] * s[0]) * idet;
	T[2][3] = (-M[2][0] * s[4] + M[2][1] * s[2] - M[2][3] * s[0]) * idet;

	T[3][0] = (-M[1][0] * c[3] + M[1][1] * c[1] - M[1][2] * c[0]) * idet;
	T[3][1] = ( M[0][0] * c[3] - M[0][1] * c[1] + M[0][2] * c[0]) * idet;
	T[3][2] = (-M[3][0] * s[3] + M[3][1] * s[1] - M[3][2] * s[0]) * idet;
	T[3][3] = ( M[2][0] * s[3] - M[2][1] * s[1] + M[2][2] * s[0]) * idet;
}
static inline void dmat4x4_orthonormalize(dmat4x4 R, dmat4x4 M)
{
	dmat4x4_dup(R, M);
	double s = 1.;
	dvec3 h;

	dvec3_norm(R[2], R[2]);
	
	s = dvec3_mul_inner(R[1], R[2]);
	dvec3_scale(h, R[2], s);
	dvec3_sub(R[1], R[1], h);
	dvec3_norm(R[2], R[2]);

	s = dvec3_mul_inner(R[1], R[2]);
	dvec3_scale(h, R[2], s);
	dvec3_sub(R[1], R[1], h);
	dvec3_norm(R[1], R[1]);

	s = dvec3_mul_inner(R[0], R[1]);
	dvec3_scale(h, R[1], s);
	dvec3_sub(R[0], R[0], h);
	dvec3_norm(R[0], R[0]);
}

static inline void dmat4x4_frustum(dmat4x4 M, double l, double r, double b, double t, double n, double f)
{
	M[0][0] = 2.f*n/(r-l);
	M[0][1] = M[0][2] = M[0][3] = 0.f;
	
	M[1][1] = 2.*n/(t-b);
	M[1][0] = M[1][2] = M[1][3] = 0.f;

	M[2][0] = (r+l)/(r-l);
	M[2][1] = (t+b)/(t-b);
	M[2][2] = -(f+n)/(f-n);
	M[2][3] = -1.f;
	
	M[3][2] = -2.f*(f*n)/(f-n);
	M[3][0] = M[3][1] = M[3][3] = 0.f;
}
static inline void dmat4x4_ortho(dmat4x4 M, double l, double r, double b, double t, double n, double f)
{
	M[0][0] = 2.f/(r-l);
	M[0][1] = M[0][2] = M[0][3] = 0.f;

	M[1][1] = 2.f/(t-b);
	M[1][0] = M[1][2] = M[1][3] = 0.f;

	M[2][2] = -2.f/(f-n);
	M[2][0] = M[2][1] = M[2][3] = 0.f;
	
	M[3][0] = -(r+l)/(r-l);
	M[3][1] = -(t+b)/(t-b);
	M[3][2] = -(f+n)/(f-n);
	M[3][3] = 1.f;
}
static inline void dmat4x4_perspective(dmat4x4 m, double y_fov, double aspect, double n, double f)
{
	/* NOTE: Degrees are an unhandy unit to work with.
	 * lindmath.h uses radians for everything! */
	double const a = 1.f / tan(y_fov / 2.f);

	m[0][0] = a / aspect;
	m[0][1] = 0.f;
	m[0][2] = 0.f;
	m[0][3] = 0.f;

	m[1][0] = 0.f;
	m[1][1] = a;
	m[1][2] = 0.f;
	m[1][3] = 0.f;

	m[2][0] = 0.f;
	m[2][1] = 0.f;
	m[2][2] = -((f + n) / (f - n));
	m[2][3] = -1.f;

	m[3][0] = 0.f;
	m[3][1] = 0.f;
	m[3][2] = -((2.f * f * n) / (f - n));
	m[3][3] = 0.f;
}
static inline void dmat4x4_look_at(dmat4x4 m, dvec3 eye, dvec3 center, dvec3 up)
{
	/* Adapted from Android's OpenGL Matrix.java.                        */
	/* See the OpenGL GLUT documentation for gluLookAt for a description */
	/* of the algorithm. We implement it in a straightforward way:       */

	/* TODO: The negation of of can be spared by swapping the order of
	 *       operands in the following cross products in the right way. */
	dvec3 f;
	dvec3_sub(f, center, eye);	
	dvec3_norm(f, f);	
	
	dvec3 s;
	dvec3_mul_cross(s, f, up);
	dvec3_norm(s, s);

	dvec3 t;
	dvec3_mul_cross(t, s, f);

	m[0][0] =  s[0];
	m[0][1] =  t[0];
	m[0][2] = -f[0];
	m[0][3] =   0.f;

	m[1][0] =  s[1];
	m[1][1] =  t[1];
	m[1][2] = -f[1];
	m[1][3] =   0.f;

	m[2][0] =  s[2];
	m[2][1] =  t[2];
	m[2][2] = -f[2];
	m[2][3] =   0.f;

	m[3][0] =  0.f;
	m[3][1] =  0.f;
	m[3][2] =  0.f;
	m[3][3] =  1.f;

	dmat4x4_translate_in_place(m, -eye[0], -eye[1], -eye[2]);
}

typedef double dquat[4];
static inline void dquat_identity(dquat q)
{
	q[0] = q[1] = q[2] = 0.f;
	q[3] = 1.f;
}
static inline void dquat_add(dquat r, dquat a, dquat b)
{
	int i;
	for(i=0; i<4; ++i)
		r[i] = a[i] + b[i];
}
static inline void dquat_sub(dquat r, dquat a, dquat b)
{
	int i;
	for(i=0; i<4; ++i)
		r[i] = a[i] - b[i];
}
static inline void dquat_mul(dquat r, dquat p, dquat q)
{
	dvec3 w;
	dvec3_mul_cross(r, p, q);
	dvec3_scale(w, p, q[3]);
	dvec3_add(r, r, w);
	dvec3_scale(w, q, p[3]);
	dvec3_add(r, r, w);
	r[3] = p[3]*q[3] - dvec3_mul_inner(p, q);
}
static inline void dquat_scale(dquat r, dquat v, double s)
{
	int i;
	for(i=0; i<4; ++i)
		r[i] = v[i] * s;
}
static inline double dquat_inner_product(dquat a, dquat b)
{
	double p = 0.f;
	int i;
	for(i=0; i<4; ++i)
		p += b[i]*a[i];
	return p;
}
static inline void dquat_conj(dquat r, dquat q)
{
	int i;
	for(i=0; i<3; ++i)
		r[i] = -q[i];
	r[3] = q[3];
}
static inline void dquat_rotate(dquat r, double angle, dvec3 axis) {
	dvec3 v;
	dvec3_scale(v, axis, sinf(angle / 2));
	int i;
	for(i=0; i<3; ++i)
		r[i] = v[i];
	r[3] = cosf(angle / 2);
}
#define dquat_norm dvec4_norm
static inline void dquat_mul_dvec3(dvec3 r, dquat q, dvec3 v)
{
/*
 * Method by Fabian 'ryg' Giessen (of Farbrausch)
t = 2 * cross(q.xyz, v)
v' = v + q.w * t + cross(q.xyz, t)
 */
	dvec3 t;
	dvec3 q_xyz = {q[0], q[1], q[2]};
	dvec3 u = {q[0], q[1], q[2]};

	dvec3_mul_cross(t, q_xyz, v);
	dvec3_scale(t, t, 2);

	dvec3_mul_cross(u, q_xyz, t);
	dvec3_scale(t, t, q[3]);

	dvec3_add(r, v, t);
	dvec3_add(r, r, u);
}
static inline void dmat4x4_from_dquat(dmat4x4 M, dquat q)
{
	double a = q[3];
	double b = q[0];
	double c = q[1];
	double d = q[2];
	double a2 = a*a;
	double b2 = b*b;
	double c2 = c*c;
	double d2 = d*d;
	
	M[0][0] = a2 + b2 - c2 - d2;
	M[0][1] = 2.f*(b*c + a*d);
	M[0][2] = 2.f*(b*d - a*c);
	M[0][3] = 0.f;

	M[1][0] = 2*(b*c - a*d);
	M[1][1] = a2 - b2 + c2 - d2;
	M[1][2] = 2.f*(c*d + a*b);
	M[1][3] = 0.f;

	M[2][0] = 2.f*(b*d + a*c);
	M[2][1] = 2.f*(c*d - a*b);
	M[2][2] = a2 - b2 - c2 + d2;
	M[2][3] = 0.f;

	M[3][0] = M[3][1] = M[3][2] = 0.f;
	M[3][3] = 1.f;
}

static inline void dmat4x4o_mul_dquat(dmat4x4 R, dmat4x4 M, dquat q)
{
/*  XXX: The way this is written only works for othogonal dmatrices. */
/* TODO: Take care of non-orthogonal case. */
	dquat_mul_dvec3(R[0], q, M[0]);
	dquat_mul_dvec3(R[1], q, M[1]);
	dquat_mul_dvec3(R[2], q, M[2]);

	R[3][0] = R[3][1] = R[3][2] = 0.f;
	R[3][3] = 1.f;
}
static inline void dquat_from_dmat4x4(dquat q, dmat4x4 M)
{
	double r=0.f;
	int i;

	int perm[] = { 0, 1, 2, 0, 1 };
	int *p = perm;

	for(i = 0; i<3; i++) {
		double m = M[i][i];
		if( m < r )
			continue;
		m = r;
		p = &perm[i];
	}

	r = sqrtf(1.f + M[p[0]][p[0]] - M[p[1]][p[1]] - M[p[2]][p[2]] );

	if(r < 1e-6) {
		q[0] = 1.f;
		q[1] = q[2] = q[3] = 0.f;
		return;
	}

	q[0] = r/2.f;
	q[1] = (M[p[0]][p[1]] - M[p[1]][p[0]])/(2.f*r);
	q[2] = (M[p[2]][p[0]] - M[p[0]][p[2]])/(2.f*r);
	q[3] = (M[p[2]][p[1]] - M[p[1]][p[2]])/(2.f*r);
}

#endif
