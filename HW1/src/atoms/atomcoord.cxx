#include <cstdlib>
#include <iostream>
#include <fstream>	
#include <vector>
#include <string>
#include "atomcoord.h"
using std::string;
using std::vector;

// Constructor
AtomCoord::AtomCoord(double x_val, double y_val, double z_val){
  x = x_val;
  y = y_val;
  z = z_val;
 
  coordVec.push_back(x);
  coordVec.push_back(y);
  coordVec.push_back(z);
}

// Public Methods
vector<double> AtomCoord::toVector(){
  return coordVec;
}

double AtomCoord::getX() const{
  return x;
}


double AtomCoord::getY() const{
  return y;
}

double AtomCoord::getZ() const{
  return z;
}

void AtomCoord::print_info() {
  std::cout << "x: " << x << ", y: " << y << ", z: " << z << std::endl;
}

std::ostream& operator<<(std::ostream& os, const AtomCoord& ac){
  os << "(" << ac.getX() << ", " << ac.getY() << ", " << ac.getZ() << ")" << "\n";
  return os;
}



