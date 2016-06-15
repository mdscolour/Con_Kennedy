#pragma once
#include "global.h"

template <class T> 
class GPoint
{
public:
	T x, y, z;
	//everything just inline b.c. short
	GPoint(){}
	GPoint(T a, T b, T c) :x(a), y(b), z(c){}
	GPoint(const GPoint& p) :x(p.x), y(p.y), z(p.z){}
	/////~point(){ cout << "deconstructed" << endl; }

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

	GPoint& operator=(const GPoint& p){ x = p.x; y = p.y; z = p.z; return(*this); }
	bool operator==(GPoint p){ return((abs(p.x - x)<1e-8 && abs(p.y - y)<1e-8 && abs(p.z - z)<1e-8)); }
	bool close(GPoint p, double preci = 1.e-5){ return((abs(p.x - x)<preci && abs(p.y - y)<preci && abs(p.z - z)<preci)); }
	T dot(GPoint p){ return(x*p.x+y*p.y+z*p.z); }
	void zero()	{x = 0;y = 0;z = 0;}
	void assign(T xx, T yy, T zz) { x = xx; y = yy; z = zz; }

	T distance(){return(T(sqrt(x*x + y*y + z*z)));}
	
	virtual T WellSeparate(){ T t = abs(x) + abs(y) + abs(z); return((t>1.e-5) ? t : 0); }
	virtual void print(FILE *fptr){ fprintf(fptr, "%f %f %f\n", x, y, z); }
	virtual void scan(FILE *fptr){ fscanf(fptr, "%f %f %f\n", &x, &y, &z); }
// 	virtual void print(FILE *fptr) = 0;
// 	virtual void scan(FILE *fptr) = 0;

	GPoint& mirror(int isym) // isym 0-7
	{
		switch (isym) {
		case 0:                        break;   // +++
		case 1:                 z = -z;  break;   // ++-
		case 2:          y = -y;         break;   // +-+
		case 3:          y = -y;  z = -z;  break;   // +--
		case 4:   x = -x;                break;   // -++
		case 5:   x = -x;         z = -z;  break;   // -+-
		case 6:   x = -x;  y = -y;         break;   // --+
		case 7:   x = -x;  y = -y;  z = -z;  break;   // ---
		default: printf("wrong mirror index, not 0-7 \n");break;
		} // end switch on isym_sign
		return (*this);
	}
	virtual void euclidean_op(GPoint* q, GPoint* p, int isym)
	{
		switch (isym % 6) {
		case 0:   x = (*q).x;     y = (*q).y;     z = (*q).z;  break;   // identity
		case 1:   x = (*q).y;     y = (*q).x;     z = (*q).z;  break;   // (xy)
		case 2:   x = (*q).x;     y = (*q).z;     z = (*q).y;  break;   // (yz)
		case 3:   x = (*q).z;     y = (*q).y;     z = (*q).x;  break;   // (xz)
		case 4:   x = (*q).z;     y = (*q).x;     z = (*q).y;  break;   // (xyz)
		case 5:   x = (*q).y;     y = (*q).z;     z = (*q).x;  break;   // (xzy) 
		default:
			printf("bad case cubic isym_perm in euclidean_op \n");
			exit(0);
			break;
		} // end switch on isym_perm
		switch (isym / 6) {
		case 0:                        break;   // +++
		case 1:                 z = -z;  break;   // ++-
		case 2:          y = -y;         break;   // +-+
		case 3:          y = -y;  z = -z;  break;   // +--
		case 4:   x = -x;                break;   // -++
		case 5:   x = -x;         z = -z;  break;   // -+-
		case 6:   x = -x;  y = -y;         break;   // --+
		case 7:   x = -x;  y = -y;  z = -z;  break;   // ---
		default:
			printf("bad case cubic isym_sign in euclidean_op \n");
			exit(0);
			break;
		} // end switch on isym_sign
		x += (*p).x;
		y += (*p).y;
		z += (*p).z;
	} // end point::euclidean_op()

	virtual GPoint& do_op(GPoint* p, int isym)
	{
		GPoint<T>* q = new GPoint<T>(x - (p->x), y - (p->y), z - (p->z));
		switch (isym % 6) {
		case 0:   x = (*q).x;     y = (*q).y;     z = (*q).z;  break;   // identity
		case 1:   x = (*q).y;     y = (*q).x;     z = (*q).z;  break;   // (xy)
		case 2:   x = (*q).x;     y = (*q).z;     z = (*q).y;  break;   // (yz)
		case 3:   x = (*q).z;     y = (*q).y;     z = (*q).x;  break;   // (xz)
		case 4:   x = (*q).z;     y = (*q).x;     z = (*q).y;  break;   // (xyz)
		case 5:   x = (*q).y;     y = (*q).z;     z = (*q).x;  break;   // (xzy) 
		default:
			printf("bad case cubic isym_perm in euclidean_op \n");
			exit(0);
			break;
		} // end switch on isym_perm
		switch (isym / 6) {
		case 0:                        break;   // +++
		case 1:                 z = -z;  break;   // ++-
		case 2:          y = -y;         break;   // +-+
		case 3:          y = -y;  z = -z;  break;   // +--
		case 4:   x = -x;                break;   // -++
		case 5:   x = -x;         z = -z;  break;   // -+-
		case 6:   x = -x;  y = -y;         break;   // --+
		case 7:   x = -x;  y = -y;  z = -z;  break;   // ---
		default:
			printf("bad case cubic isym_sign in euclidean_op \n");
			exit(0);
			break;
		} // end switch on isym_sign
		x += (*p).x;
		y += (*p).y;
		z += (*p).z;
		delete q;
		return(*this);
	} // end point::do_op()

	virtual GPoint get_op(GPoint* p, int isym)
	{
		GPoint<T>* q = new GPoint<T>(x - (p->x), y - (p->y), z - (p->z));
		T x, y, z;
		switch (isym % 6) {
		case 0:   x = (*q).x;     y = (*q).y;     z = (*q).z;  break;   // identity
		case 1:   x = (*q).y;     y = (*q).x;     z = (*q).z;  break;   // (xy)
		case 2:   x = (*q).x;     y = (*q).z;     z = (*q).y;  break;   // (yz)
		case 3:   x = (*q).z;     y = (*q).y;     z = (*q).x;  break;   // (xz)
		case 4:   x = (*q).z;     y = (*q).x;     z = (*q).y;  break;   // (xyz)
		case 5:   x = (*q).y;     y = (*q).z;     z = (*q).x;  break;   // (xzy) 
		default:
			printf("bad case cubic isym_perm in euclidean_op \n");
			exit(0);
			break;
		} // end switch on isym_perm
		switch (isym / 6) {
		case 0:                        break;   // +++
		case 1:                 z = -z;  break;   // ++-
		case 2:          y = -y;         break;   // +-+
		case 3:          y = -y;  z = -z;  break;   // +--
		case 4:   x = -x;                break;   // -++
		case 5:   x = -x;         z = -z;  break;   // -+-
		case 6:   x = -x;  y = -y;         break;   // --+
		case 7:   x = -x;  y = -y;  z = -z;  break;   // ---
		default:
			printf("bad case cubic isym_sign in euclidean_op \n");
			exit(0);
			break;
		} // end switch on isym_sign
		x += (*p).x;
		y += (*p).y;
		z += (*p).z;
		delete q;
		return(GPoint(x,y,z));
	} // end point::do_op()
};

// class fPoint :public GPoint<float>
// {
// public:
// 	fPoint(){}
// 	fPoint(float a, float b, float c) :GPoint(a, b, c){}
// 	fPoint(const GPoint& p):GPoint(p.x,p.y,p.z){}
// 	virtual void print(FILE *fptr){ fprintf(fptr, "%lf %lf %lf\n", x, y, z); }
// };
// class iPoint :public GPoint<int>
// {
// public:
// 	iPoint(){}
// 	iPoint(int a, int b, int c) :GPoint(a, b, c){}
// 	iPoint(const GPoint& p) :GPoint(p.x, p.y, p.z){}
// 	virtual void print(FILE *fptr){ fprintf(fptr, "%d %d %d\n", x, y, z); }
// };

class Matrix
{
public:
	GPoint<double> row1;
	GPoint<double> row2;
	GPoint<double> row3;
	Matrix(){}
	Matrix(char* key){ if (key == "identity")this->identity(); }
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

	void print(FILE *fptr){ row1.print(fptr); row2.print(fptr); row3.print(fptr); }
	Matrix& identity(){ row1 = GPoint<double>(1, 0, 0); row2 = GPoint<double>(0, 1, 0); row3 = GPoint<double>(0, 0, 1); return(*this); }
	GPoint<double> dot(GPoint<double> v){ return GPoint<double>(row1.dot(v), row2.dot(v), row3.dot(v)); }
	Matrix dot(Matrix x){ return Matrix(row1.dot(x.col1()), row1.dot(x.col2()), row1.dot(x.col3()), row2.dot(x.col1()), row2.dot(x.col2()), row2.dot(x.col3()), row3.dot(x.col1()), row3.dot(x.col2()), row3.dot(x.col3())); }
	//Sphere dot(Sphere v){return GPoint<double>(row1.dot(v),row2.dot(v),row3.dot(v));}

};
class RMatrix : public Matrix
{
public:
	RMatrix(){}
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
		double theta1 = (isym % RADIANNUM) * 2.0 * M_PI / RADIANNUM; isym /= RADIANNUM;
		double theta2 = ((isym % RADIANNUM)+1) * 2.0 * M_PI / (RADIANNUM+1); isym /= RADIANNUM;
		double theta3 = ((isym % RADIANNUM)+1) * 2.0 * M_PI / (RADIANNUM+1); isym /= RADIANNUM;
		new(this) RMatrix(isym % ROTNUM, theta1, theta2, theta3);
	}
	//RMatrix Rotbyind(long isym, double theta){ return(RMatrix(isym, theta1, theta2)); }
	Matrix rotxy(double theta){ return Matrix(cos(theta), -sin(theta), 0, sin(theta), cos(theta), 0, 0, 0, 1); }
	Matrix rotxz(double theta){ return Matrix(cos(theta), 0, -sin(theta), 0, 1, 0, sin(theta), 0, cos(theta)); }
	Matrix rotyz(double theta){ return Matrix(1, 0, 0, 0, cos(theta), -sin(theta), 0, sin(theta), cos(theta)); }
};

class RefMatrix : public Matrix
{
public:
	RefMatrix(){}
	RefMatrix(Matrix x) :Matrix(x.row1, x.row2, x.row3){}
	RefMatrix(int isym)
	{
		double x=1, y=1, z=1;
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
		new(this) Matrix(x,0,0,0,y,0,0,0,z);
	}
};

class OpMatrix :public Matrix
{
public:
	OpMatrix(){}
	OpMatrix(Matrix x) :Matrix(x.row1, x.row2, x.row3){}
	OpMatrix(int isym)
	{
		RMatrix temp(isym);
		isym /= RADIANNUM;
		isym /= RADIANNUM;
		isym /= RADIANNUM;
		new(this) OpMatrix(RefMatrix(isym / ROTNUM).dot(temp));
	}
};


class Sphere :public GPoint<double>
{
public:
	double r;
	Sphere() :r(RADIUS){}
	Sphere(double a, double b, double c) :GPoint(a, b, c), r(RADIUS){}
	Sphere(double a, double b, double c, double d) :GPoint(a, b, c), r(d){}
	Sphere(const GPoint& p) :GPoint(p.x, p.y, p.z),r(RADIUS){}
	virtual void print(FILE *fptr){ fprintf(fptr, "%lf %lf %lf %lf \n", x, y, z, r); }
	virtual void scan(FILE *fptr){ fscanf(fptr, "%lf %lf %lf %lf \n", &x, &y, &z, &r); }
	virtual double WellSeparate(){ double t = sqrt(x*x + y*y + z*z); return((t > 2*r) ? t : 0); }
	virtual Sphere get_op(GPoint<double>* ref, OpMatrix* op)
	{
		return( op->dot(*this)+(*ref));
	}
	virtual void euclidean_op(GPoint<double>* p,GPoint<double>* ref, OpMatrix* op)
	{
		(*this) = op->dot(*p) + (*ref);
	}
};



