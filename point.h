/**************************************************************************
The basic of the monomer of the model. Containing this class:
GPoint<T> --- 3d point and general operation

Matrix --- 3d matrix and general operation

RMatrix --- rotation matrix, RMatrix() will return a random one.

RefMatrix --- reflection matrix, RefMatrix() will return a random one.

This two are more important:

OpMatrix --- matrix of opeation, currently is just combination of rotation and reflection.
For safety, OpMatrix() will return a blank matrix. 
OpMatrix("rand") will return a random one, where most of the case we want.

Sphere --- hard shpere.

In class Walk, (a-b).WellSeparate() will be called to judge whether a and b overlapped, and how to increase the index.
When the return is 0, they are overlapped and deny the pivot. 
When not 0 the index will be increased by the amount of the return, i.e. within (return-1) steps will be regarded as safe.

In class Walk, a.euclidean_op(&b,&c,&op) will be called to realize the function: op.dot(b)+c and store in a.

Override or hide the above two function if a new model is made. (no vitural function is needed in this version.)
**************************************************************************/

#pragma once
#include "global.h"

class Sphere;
class OpMatrix;
// Template class of of point, all ordinary operation of 3x3 point is available
template <class T> 
class GPoint
{
public:
	T x, y, z;
	//everything just inline
	GPoint(){}
	GPoint(T a, T b, T c) :x(a), y(b), z(c){}
	//GPoint(const Sphere& p) :x(p.x), y(p.y), z(p.z){printf("aaa done once\n");}
	//GPoint(const GPoint& p) :x(p.x), y(p.y), z(p.z){}
	GPoint& operator=(const GPoint& p){ x = p.x; y = p.y; z = p.z; return(*this); }

	T coord_x(){ return(x); }
	T coord_y(){ return(y); }
	T coord_z(){ return(z); }
	GPoint operator+(GPoint p){ return GPoint(x + p.x, y + p.y, z + p.z); }
	GPoint operator-(GPoint p){ return GPoint(x - p.x, y - p.y, z - p.z); }
	GPoint operator*(T c){ return GPoint(x*c, y*c, z*c); }
	GPoint operator/(T c){ return GPoint(x/c, y/c, z/c); }
	GPoint& operator+=(GPoint p){ x += p.x; y += p.y; z += p.z; return *this; }
	GPoint& operator-=(GPoint p){ x -= p.x; y -= p.y; z -= p.z; return *this; }
	GPoint& operator*=(T c){ x *= c; y *= c; z *= c; return *this; }
	GPoint& operator/=(T c){ x /= c; y /= c; z /= c; return *this; }

	bool operator==(GPoint p){ return((abs(p.x - x)<1e-8 && abs(p.y - y)<1e-8 && abs(p.z - z)<1e-8)); }
	bool close(GPoint p, double preci = 1.e-5){ return((abs(p.x - x)<preci && abs(p.y - y)<preci && abs(p.z - z)<preci)); }

	T dot(GPoint p){ return(x*p.x+y*p.y+z*p.z); }
	inline T dot(Sphere p);
	
	void zero()	{x = 0;y = 0;z = 0;}
	//void assign(T xx, T yy, T zz) { x = xx; y = yy; z = zz; }

	T norm(){return(T(sqrt(x*x + y*y + z*z)));}
	inline GPoint rotation(OpMatrix* op);
};

// The real model which is a hard sphere
class Sphere
{
public:
	double x,y,z,r,k;
	Sphere(){}
	Sphere(double a, double b, double c, double d, double e) :x(a),y(b),z(c),r(d),k(e){}
	//Sphere(const GPoint& p) :GPoint(p.x, p.y, p.z),r(RADIUS){}

	Sphere operator+(GPoint<double> p){ return Sphere(x + p.x, y + p.y, z + p.z, r, k); }
	//Sphere operator-(GPoint<double> p){ return Sphere(x - p.x, y - p.y, z - p.z, r); }
	GPoint<double> operator-(Sphere p){ return GPoint<double>(x - p.x, y - p.y, z - p.z); }
	GPoint<double> topoint(){return GPoint<double>(x,y,z);}
	
	void assign(double xx, double yy, double zz, double rr,double kk) { x = xx; y = yy; z = zz; r = rr; k = kk; }
	void print(FILE *fptr){ fprintf(fptr, "%lf %lf %lf %lf %lf \n", x, y, z, r, k); }
	void scan(FILE *fptr){ fscanf(fptr, "%lf %lf %lf %lf %lf \n", &x, &y, &z, &r, &k); }
	int WellSeparate(const Sphere& s){ double t = sqrt((x-s.x)*(x-s.x) + (y-s.y)*(y-s.y) + (z-s.z)*(z-s.z)); return((t > s.r+this->r) ? int(t+(1-s.r+this->r)) : 0); }
	double distance(const Sphere& s){return(sqrt((x-s.x)*(x-s.x) + (y-s.y)*(y-s.y) + (z-s.z)*(z-s.z)));}
	//double norm(){return(sqrt(x*x + y*y + z*z));}
	inline void euclidean_op(Sphere* p, GPoint<double>* ref, OpMatrix* op);
};

// 3x3 matrix
class Matrix
{
public:
	GPoint<double> row1;
	GPoint<double> row2;
	GPoint<double> row3;
	Matrix(){}
	Matrix(const char* key){ if (key == "identity")this->identity(); }
	Matrix(GPoint<double> a, GPoint<double> b, GPoint<double> c) :row1(a), row2(b), row3(c){}
	Matrix(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3) :row1(x1, y1, z1), row2(x2, y2, z2), row3(x3, y3, z3){}
	//Matrix(RMatrix x) :row1(x.row1), row2(x.row2), row3(x.row3){}

	GPoint<double> col1(){ return GPoint<double>(row1.coord_x(), row2.coord_x(), row3.coord_x()); }
	GPoint<double> col2(){ return GPoint<double>(row1.coord_y(), row2.coord_y(), row3.coord_y()); }
	GPoint<double> col3(){ return GPoint<double>(row1.coord_z(), row2.coord_z(), row3.coord_z()); }
	Matrix& operator=(const Matrix& x){ row1 = x.row1;  row2 = x.row2; row3 = x.row3; return(*this); }
	bool operator==(Matrix m){ return((row1 == m.row1 && row2 == m.row2 && row3 == m.row3)); }
	bool close(Matrix m, double preci = 0.01){ return((row1.close(m.row1, preci) && row2.close(m.row2, preci) && row3.close(m.row3, preci))); }
	Matrix inv();

	//void print(FILE *fptr){ row1.print(fptr); row2.print(fptr); row3.print(fptr); }
	Matrix& identity(){ row1 = GPoint<double>(1, 0, 0); row2 = GPoint<double>(0, 1, 0); row3 = GPoint<double>(0, 0, 1); return(*this); }
	
	GPoint<double> dot(GPoint<double> v){ return GPoint<double>(row1.dot(v), row2.dot(v), row3.dot(v)); }	
	Sphere dot(Sphere v){ return Sphere(row1.dot(v), row2.dot(v), row3.dot(v), v.r, v.k); }
	Matrix dot(Matrix x){ return Matrix(row1.dot(x.col1()), row1.dot(x.col2()), row1.dot(x.col3()), row2.dot(x.col1()), row2.dot(x.col2()), row2.dot(x.col3()), row3.dot(x.col1()), row3.dot(x.col2()), row3.dot(x.col3())); }
	//Sphere dot(Sphere v){return GPoint<double>(row1.dot(v),row2.dot(v),row3.dot(v));}

};

// 3x3 rotation matrix
class RMatrix : public Matrix
{
public:
	RMatrix()
	{
		int plane = int(RNG_NAME() * 3);
		double theta = 2 * M_PI*RNG_NAME();
		if (plane==0)new(this) Matrix(cos(theta), -sin(theta), 0, sin(theta), cos(theta), 0, 0, 0, 1);
		else if (plane == 1)new(this) Matrix(cos(theta), 0, -sin(theta), 0, 1, 0, sin(theta), 0, cos(theta));
		else if (plane == 2)new(this) Matrix(1, 0, 0, 0, cos(theta), -sin(theta), 0, sin(theta), cos(theta));
		else printf("Error in constructing RMatrix\n");
	}
	//RMatrix(){ new(this)RMatrix(int(15 * 100*100*100 * RNG_NAME())); }
	RMatrix(Matrix x) :Matrix(x.row1, x.row2, x.row3){}

	// rotation index 
	// 0 xy
	// 1 xy + xz
	// 2 xy + yz
	// 3 yz
	// 4 yz + xy
	// 5 yz + xz
	// 6 xz
	// 7 xz + xy
	// 8 xz + yz
	// 9 xy + xz + yz:123
	// 10 xy yz xz  132,
	//11 xz xy yz   213, 
	//12 yz xz xy	231, 
	//13 xz xy yz   312, 
	//14 xz yz xy	321.
	RMatrix(int isym, double theta1, double theta2, double theta3)
	{
		Matrix t1("identity"), t2("identity"), t3("identity");
		switch (isym)
		{
		case 0:{ t1 = rotxy(theta1);break; }
		case 1:{ t1 = rotxy(theta1);t2 = rotxz(theta2) ;break; }
		case 2:{ t1 = rotxy(theta1); t2 = rotyz(theta2);break; }
		case 3:{ t1 = rotyz(theta1);break; }
		case 4:{ t1 = rotyz(theta1); t2 = rotxy(theta2); break; }
		case 5:{ t1 = rotyz(theta1); t2 = rotxz(theta2); break; }
		case 6:{ t1 = rotxz(theta1); break; }
		case 7:{ t1 = rotxz(theta1); t2 = rotxy(theta2); break; }
		case 8:{ t1 = rotxz(theta1); t2 = rotyz(theta2); break; }
		case 9:{ t1 = rotxy(theta1); t2 = rotxz(theta2); t3 = rotyz(theta3); break; }
		case 10:{ t1 = rotxy(theta1); t2 = rotyz(theta2); t3 = rotxz(theta3); break; }
		case 11:{ t1 = rotxz(theta1); t2 = rotxy(theta2); t3 = rotyz(theta3); break; }
		case 12:{ t1 = rotyz(theta1); t2 = rotxz(theta2); t3 = rotxy(theta3); break; }
		case 13:{ t1 = rotxz(theta1); t2 = rotxy(theta2); t3 = rotyz(theta3); break; }
		case 14:{ t1 = rotxz(theta1); t2 = rotyz(theta2); t3 = rotxy(theta3); break; }
		default: printf("bad case in line_initialize() d\n"); exit(0); break;
		}
		new(this) RMatrix(t1.dot(t2.dot(t3)));
	}
	RMatrix(int isym)
	{
		double theta1 = (isym % 100) * 2.0 * M_PI / 100; isym /= 100;
		double theta2 = ((isym % 100)+1) * 2.0 * M_PI / (100+1); isym /= 100;
		double theta3 = ((isym % 100)+1) * 2.0 * M_PI / (100+1); isym /= 100;
		new(this) RMatrix(isym, theta1, theta2, theta3);
	}
	//RMatrix Rotbyind(long isym, double theta){ return(RMatrix(isym, theta1, theta2)); }
	Matrix rotxy(double theta){ return Matrix(cos(theta), -sin(theta), 0, sin(theta), cos(theta), 0, 0, 0, 1); }
	Matrix rotxz(double theta){ return Matrix(cos(theta), 0, -sin(theta), 0, 1, 0, sin(theta), 0, cos(theta)); }
	Matrix rotyz(double theta){ return Matrix(1, 0, 0, 0, cos(theta), -sin(theta), 0, sin(theta), cos(theta)); }
};

// 3x3 reflection matrix
class RefMatrix : public Matrix
{
public:
	RefMatrix()
	{
		int isym = int(8 * RNG_NAME());
		double x = 1, y = 1, z = 1;
		switch (isym)
		{
		case 0:                        break;   // +++
		case 1:                 z = -1;  break;   // ++-
		case 2:          y = -1;         break;   // +-+
		case 3:          y = -1;  z = -1;  break;   // +--
		case 4:   x = -1;                break;   // -++
		case 5:   x = -1;         z = -1;  break;   // -+-
		case 6:   x = -1;  y = -1;         break;   // --+
		case 7:   x = -1;  y = -1;  z = -1;  break;   // ---
		default: printf("wrong mirror index, not 0-7 \n"); break;
		} // end switch on isym_sign
		new(this) Matrix(x, 0, 0, 0, y, 0, 0, 0, z);
	}
	RefMatrix(Matrix x) :Matrix(x.row1, x.row2, x.row3){}
	RefMatrix(int isym)
	{
		double x = 1, y = 1, z = 1;
		switch (isym)
		{
		case 0:                        break;   // +++
		case 1:                 z = -1;  break;   // ++-
		case 2:          y = -1;         break;   // +-+
		case 3:          y = -1;  z = -1;  break;   // +--
		case 4:   x = -1;                break;   // -++
		case 5:   x = -1;         z = -1;  break;   // -+-
		case 6:   x = -1;  y = -1;         break;   // --+
		case 7:   x = -1;  y = -1;  z = -1;  break;   // ---
		default: printf("wrong mirror index, not 0-7 \n"); break;
		} // end switch on isym_sign
		new(this) Matrix(x, 0, 0, 0, y, 0, 0, 0, z);
	}
};

// The real used operation matrix, which is nothing but reflection matrix dot rotation matrix.
// can only be constructed by matrix of a integer, the integer must be 0 to (num_rotation*num_reflection-1)*RADIANNUM*RADIANNUM*RADIANNUM
class OpMatrix :public Matrix
{
public:
	OpMatrix(){}
	OpMatrix(const char* key){ if (key == "rand")new(this) OpMatrix(RefMatrix().dot(RMatrix())); }
	OpMatrix(Matrix x) :Matrix(x.row1, x.row2, x.row3){}
	OpMatrix(int isym)
	{
		RMatrix temp(isym);
		isym /= 100;
		isym /= 100;
		isym /= 100;
		new(this) OpMatrix(RefMatrix(isym / 15).dot(temp));
	}
};

template <class T> 
inline T GPoint<T>::dot(Sphere p){ return(x*p.x+y*p.y+z*p.z); }

template <class T> 
inline GPoint<T> GPoint<T>::rotation(OpMatrix* op)
{
	return(op->dot(*this));
}

inline void Sphere::euclidean_op(Sphere* p, GPoint<double>* ref, OpMatrix* op)
{
	(*this) = op->dot(*p) + (*ref);
}

