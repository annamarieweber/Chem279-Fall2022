/*
  UC Berkeley - MSSE Program
  Chem 279 
  Fall 2022
  This file, atomcoord.h, contains an API to the AtomCoord class.
  being used in Chem 279.
*/
#ifndef ATOM_COORD
#define ATOM_COORD
#include <cstdlib>
#include <iostream>
#include <fstream>	
#include <vector>
#include <string>

using std::string;
using std::vector;

/**
 * AtomCoord - class for holding the coordinate an atom is located at
 *
 **/
class AtomCoord{
  private:
    double x;
    double y;
    double z;
    vector<double> coordVec;

  public:
    /**
     * @brief AtomCoord constructor
     *
     * @param x_val double representing the x atomic coordinate
     * @param y_val double representing the y atomic coordinate
     * @param z_val double representing the z atomic coordinate
     *
     **/
    AtomCoord(double x_val, double y_val, double z_val);
    
    
    /** 
     * @brief converts members of coord to vector of double 
     *
     * @return vector<double> - a vector of lenght 3 in the form [x,y,z]
     *
     * \todo {in a perfect world I likely want to overload amadillo ops so that you can do matrix operations on vectors of(AtomCoord)}
     **/
    vector<double> toVector();

    /**
     * @brief gets the x coordinate of the vector
     *
     * @return double - x coordinate
     **/
    double getX() const;


    /**
     * @brief gets the y coordinate of the vector
     *
     * @return double - y coordinate
     **/
    double getY() const;

    /**
     * @brief gets the z coordinate of the vector
     *
     * @return double - z coordinate
     **/
    double getZ() const;

    /**
     * @brief print the coordintates info
     **/
    void print_info();


    /**
     * <<
     *
     * overloaded outputstream operator for AtomCoord class
     *
     * @param os ostream - the output stream
     * @param a AtomCoord- a reference to the Object being printed
     **/
    friend std::ostream& operator<<(std::ostream& os, const AtomCoord& ac);
};
#endif

