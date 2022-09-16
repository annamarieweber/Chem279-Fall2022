#include <atomic>
#include <cstdlib>
#include <iostream>
#include <fstream>	
#include <vector>
#include <string>
#include <map>
#include <cstdlib>
#include "atomcoord.h"
#include "atominfos.h"
#include "constants.h"
#include "clustermap.h"
#include "baseexception.h"

using std::string;
using std::vector;
using constants::ATOM_SYMBOLS;

/**
 * ClusterMap - class for holding the information about All atoms in a cluster
 *
 **/
ClusterMap::ClusterMap(){
}

AtomInfos ClusterMap::getAtomInfos(int num) const{
  try{
    string sym = ATOM_SYMBOLS[num-1];
    if(atomTypes.count(sym) < 1 ){
      std::cout << "Error: Cluster has no Atoms with number: " << num << std::endl;
    }
    else{
      return atomTypes.at(sym);
    }
  }catch(...){
    std::cout << "somethign went wrong" << std::endl;
  }
  return AtomInfos(0);
  
}


void ClusterMap::addCoord(int num, double x_val, double y_val, double z_val){
  try{
    string sym = ATOM_SYMBOLS[num-1];
    //todo: remove hard coded check for gold
    if(sym != "Au"){
      throw BaseException("DisallowedAtomType: Only Au atoms are allowed.");
    }
    else{
      if(atomTypes.count(sym) < 1){
	atomTypes.insert(std::pair<string,AtomInfos>(sym,AtomInfos(num)));
      }
      atomTypes.at(sym).addCoord(x_val, y_val, z_val);
    }
  }
  catch (BaseException e){
    e.displayError();
  }
  catch (...){
    std::cout << "Something went wrong" << std::endl;
  }
}

std::ostream& operator<<(std::ostream& os, const ClusterMap& c){

   for(const auto& atomPair : c.atomTypes){
     os << "Atom Symbol: " << atomPair.first << "\n";
     os << atomPair.second;
   }

   return os;
}



