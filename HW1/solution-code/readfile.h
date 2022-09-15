#ifndef READFILE_H
#define READFILE_H

#include <iostream>
#include <vector>


struct Atom
{
  int AtomicNumber;
  double r[3];
  friend std::ostream &operator<<(std::ostream &os, const Atom p)
  {
    os << p.AtomicNumber << "(" << p.r[0] << ", "
       << p.r[1] << ", " << p.r[2] << ")";
    return os;
  }
};

void Readatomsfromfile(std::vector<Atom> &Atoms, std::string &fname);


#endif // READFILE_H