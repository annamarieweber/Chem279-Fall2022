/*
  UC Berkeley - MSSE Program
  Chem 279 
  Fall 2022
  This file, clustermap.h, contains an API to the ClusterMap class.
  being used in Chem 279.
*/
#ifndef CLUSTER_MAP
#define CLUSTER_MAP

#include <exception>
#include <string>
#include <map>
#include "constants.h"
#include "atominfos.h"

using std::string;
using std::map;
using constants::ATOM_SYMBOLS;

class ClusterMap{
  private:
    map<string, AtomInfos> atomTypes;
  public:
    /**
     * ClusterMap constructor
     **/
    ClusterMap();

    /**
     * @brief gets AtomInfos for the specified symbol if it exists
     *
     * @param num int representing the atomic number for the element to get AtomInfos for
     * @return AtomInfos
     **/
    AtomInfos getAtomInfos(int num) const;

    /**
     * @brief adds a neew atomic coordinate
     *
     * @detail adds a new atomic coordinate if the atomis not already in the ClusterMap creates a new maping and adds coordinate to that
     *
     * @param num int representing the atomic number of the atom being added
     * @param x_val double x coordinate
     * @param y_val double y coordinate
     * @param z_val doubel z coordinte
     **/
    void addCoord(int num, double x_val, double y_val, double z_val);

    friend std::ostream& operator<<(std::ostream& os, const ClusterMap& c);
};
#endif
