#pragma once
#include <cmath>


#define PI 3.14159
// Simple vector struct to do vector operations.


struct vector2D{
	vector2D(){};
	vector2D(double X, double Y): x(X), y(Y){};
	double x;
	double y;
};

struct vector3D:public vector2D {
	//vector3D(){};
	vector3D(double X=0, double Y=0, double Z=0):vector2D(X,Y),z(Z){};
	double z;
};

inline vector3D operator+(vector3D const & A,vector3D const & B  ){
	vector3D result;
	result.x = A.x+B.x;
	result.y = A.y+B.y;
	result.z = A.z+B.z;
	return result;
}

inline vector3D operator-(vector3D const & A,vector3D const & B  ){
	vector3D result;
	result.x = A.x-B.x;
	result.y = A.y-B.y;
	result.z = A.z-B.z;
	return result;
}

inline vector3D operator/(vector3D const & A,double const & B  ){
	vector3D result;
	result.x = A.x/B;
	result.y = A.y/B;
	result.z = A.z/B;
	return result;
}

inline vector3D operator*(vector3D const & A,double const & B  ){
	vector3D result;
	result.x = A.x*B;
	result.y = A.y*B;
	result.z = A.z*B;
	return result;
}

inline vector3D operator*(double const & B, vector3D const & A  ){
	vector3D result;
	result.x = A.x*B;
	result.y = A.y*B;
	result.z = A.z*B;
	return result;
}

// Vector 2D
 inline vector3D operator+(vector2D const & A,vector2D const & B  ){
	vector3D result;
	result.x = A.x+B.x;
	result.y = A.y+B.y;
	return result;
}

inline vector2D operator-(vector2D const & A,vector2D const & B  ){
	vector2D result;
	result.x = A.x-B.x;
	result.y = A.y-B.y;
	return result;
}

inline vector2D operator*(vector2D const & A,double const & B  ){
	vector2D result;
	result.x = A.x*B;
	result.y = A.y*B;
	return result;
}

inline vector2D operator/(vector2D const & A,double const & B  ){
	vector2D result;
	result.x = A.x/B;
	result.y = A.y/B;
	return result;
}

inline vector2D operator*(double const & B, vector2D const & A  ){
	vector2D result;
	result.x = A.x*B;
	result.y = A.y*B;
	return result;
}




// cross product between two vectors
inline vector3D cross(vector3D const & A,vector3D const & B  ){
	vector3D result;
	result.x = A.y*B.z - A.z*B.y ;
	result.y = -(A.x*B.z-A.z*B.x);
	result.z = A.x*B.y-A.y*B.x;
	return result;
}

// dot product between two vectors
inline double dot(vector3D const & A,vector3D const & B  ){
	return A.x*B.x+A.y*B.y+A.z*B.z;
}

inline double findAngle(vector3D p1, vector3D p2){
	 double angle = std::acos(dot(p1, p2));
	 if (p1.x*p2.y < p1.y*p2.x)
	 {
		 angle = 2 * PI -angle;
	 }
	 return angle;
};


// Magnitud of a vector
inline double norm(vector3D const & A  ){
	return sqrt(A.x*A.x+A.y*A.y+A.z*A.z  );
}

