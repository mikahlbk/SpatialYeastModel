// coord.cpp
//*************************************************

// forward dependencies
//*************************************************

// include dependencies
#include <omp.h>
#include <iostream>
#include <unistd.h>
#include "coord.h"
#include <math.h>
//************************************************

using namespace std;

// Coord member functions

// Constructors

Coord::Coord() {
	x = 0;
	y = 0;
}

Coord::Coord(double x, double y) {
	this->x = x;
	this->y = y;
}

Coord::Coord(const Coord& c) {
	x = c.get_X();
	y = c.get_Y();
}

//Getter Functions

double Coord::get_X() const {return x;}
double Coord::get_Y() const {return y;}

//Overloading Operators

void Coord::operator=(const Coord& c){
	x = c.get_X();
	y = c.get_Y();
}
Coord Coord::operator+(const Coord& c) const{
	Coord q(x+c.get_X(),y+c.get_Y());
	return q;
}
Coord Coord::operator+(const double& d) const{
	Coord q(x+d,y+d);
	return q;
}
Coord Coord::operator-(const Coord& c) const{
	Coord q(x-c.get_X(),y-c.get_Y());
	return q;
}
Coord Coord::operator-(const double& d) const{
	Coord q(x-d,y-d);
	return q;
}
Coord Coord::operator/(const double& d) const{
	Coord q( x/d, y/d);
	return q;
}
Coord Coord::operator*(const double& d) const{
	Coord q( x*d, y*d);
	return q;
}
void Coord::operator+=(const Coord&c) {
	x += c.get_X();
	y += c.get_Y();
	return;
}
void Coord::operator+=(const double&d) {
	x += d;
	y += d;
	return;
}
void Coord::operator-=(const Coord&c) {
	x -= c.get_X();
	y -= c.get_Y();
	return;
}
void Coord::operator-=(const double&d) {
	x -= d;
	y -= d;
	return;
}
bool Coord::operator==(const Coord& c) {
	if((x == c.get_X())&&(y==c.get_Y())){
		return true;
	}
	else {
		return false;
	}
}

bool Coord::operator!=(const Coord& c) {
	if((x != c.get_X())||(y!=c.get_Y())){
		return true;
	}
	else {
		return false;
	}
}

//Higher Level Math

double Coord::length() const {
	return sqrt( (x*x) + (y*y));
}
double Coord::dot(const Coord& c) const{
	return ( (x*c.get_X()) + (y*c.get_Y()));
}
double Coord::cross(const Coord& c) const{
	return (x*c.get_Y())-(y*c.get_X());
}
// Display Functions
ostream& operator<<(ostream& os, const Coord& c) {
	os << '(' << c.get_X() << ',' << c.get_Y() << ')';
	return os;
}
