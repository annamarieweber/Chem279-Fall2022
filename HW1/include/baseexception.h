#ifndef BASE_EXCEPTION
#define BASE_EXCEPTION
#include <exception>
#include <string>
#include <iostream>

using std::string;

class BaseException : public std::exception {
  private:
    string message;
  public:
    BaseException(string msg);
    void displayError();
};

#endif
