#include <cstdlib>
#include <iostream>
#include <fstream>	
#include <vector>
#include <string>
#include "atomcoord.h"
#include "atominfos.h"
using std::string;
using std::vector;

/**
 * AtomInfos - class for holding the atomic nnumber and all atomic coordinates for a given atom type
 *
 **/
AtomInfos::AtomInfos(int atomicNumVal){
  atomicNum = atomicNumVal;
}

void AtomInfos::addCoord(double x, double y, double z){
  AtomCoord new_coord(x,y,z);
  coords.push_back(new_coord);
}

int AtomInfos::getAtomicNum() const{
 return atomicNum;
}

vector<AtomCoord> AtomInfos::getCoords() const{
 return coords;
}

std::ostream& operator<<(std::ostream& os, const AtomInfos& a){

   for(auto c : a.getCoords()){
     os << c;
   }

   return os;
}


