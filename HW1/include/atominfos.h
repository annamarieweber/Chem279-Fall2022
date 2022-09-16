#ifndef ATOM_INFOS
#define ATOM_INFOS
#include "atomcoord.h"
#include <vector>

/**
 * AtomInfos - class for holding the atomic nnumber and all atomic coordinates for a given atom type
 *
 **/
class AtomInfos{
  private:
    vector<AtomCoord> coords;
    int atomicNum;
  public:
    /**
     * AtomInfos constructor
     *
     * @param atomicNumVal int representing the atomic number of an atom
     **/
    AtomInfos(int atomicNumVal);

    /**
     * addCoord adds a new coordinte to the list of atom coordinates 
     * 
     * @param x double - the x coordinate 
     * @param y double - the y coordinate
     * @param z double - the z coordinate
     **/
     void addCoord(double x, double y, double z);

     /**
      * getAtomicNum gets the atomic number
      *
      * @return int - the atomic Number
      **/
     int getAtomicNum() const;

     /**
      * getCoords gets the vector of coordinates for the atom type
      *
      * @return vector<AtomCoord> the list of coordinates
      **/
     vector<AtomCoord> getCoords() const;

     /**
      * <<
      *
      * overloaded outputstream operator for AtomInfos class
      *
      * @param os ostream - the output stream
      * @param a AtomInfos- a reference to the Object being printed
      **/
     friend std::ostream& operator<<(std::ostream& os, const AtomInfos& a);
};
#endif
