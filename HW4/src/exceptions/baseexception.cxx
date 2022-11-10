#include "baseexception.h"
#include <string>
#include <iostream>

using std::string;

BaseException::BaseException(string msg){
  message = msg;
}

void BaseException::displayError(){
  std::cout << message << std::endl;
}

