#import <armadillo>
#import <iostream>
#import <cmath>
using arma::vec;
using std::string;

class Function
{
public:
  template <typename T>
  T operator()(T x)
  {
    return 1 * x;
  }

  template <typename T>
  T first_deriv(T x)
  {
    return 1 * x;
  }

  template <typename T>
  T second_deriv(T x)
  {
    return 1 * x;
  }
};

class Gaussian
{
private:
  double _x_a;
  double _l_a;
  double _alpha;

public:
  /**
   * copy constructor for a gaussian
   **/
  void operator = (const Gaussian &g){
    _x_a = g._x_a;
    _l_a = g._l_a;
    _alpha = g._alpha;
  }

  Gaussian(){};
  /**
   * constructor for a gaussian
   **/
  Gaussian(double x_a, double l_a, double alpha)
  {
    _x_a = x_a;
    _l_a = l_a;
    _alpha = alpha;
  }

  /**
   * evaluation of the Gaussian Function = (x − XA)
   **/
  template <typename T>
  T operator()(T x)
  {
    T res = pow(x - _x_a, _l_a) * exp(-1.0 * _alpha * pow(x - _x_a, 2));
    return res;
  };

  /**
   * first derivative of the gaussian
   * e^(-alpha (x - _x_a)^2) (x - _x_a)^(_l_a - 1) (_l_a - 2 alpha (x - _x_a)^2)
   **/
  template <typename T>
  T first_deriv(T x)
  {
    return exp(-1.0 * _alpha * pow(x - _x_a, 2)) * pow(x - _x_a, _l_a - 1) * (_l_a - 2 * _alpha * pow(x - _x_a, 2));
  }

  /**
   * second derivative of the gaussian
   * e^(-_alpha (-_x_a + x)^2) (-_x_a + x)^(_l_a - 2) (4 _alpha^2 (_x_a - x)^4 - 2 _alpha (2 _l_a + 1) (-_x_a + x)^2 + _l_a (-1 + _l_a))
   **/
  template <typename T>
  T second_deriv(T x)
  {
    return exp(-1.0 * _alpha * pow(-1.0 * (_x_a + x), 2) * pow(-1.0 * _x_a + x, _l_a - 2)) * (4 * pow(_alpha, 2) * pow(_x_a - x, 4) - 2 * _alpha * (2 * _l_a + 1) * (-1 * pow(_x_a + x, 2)) + _l_a * (-1 + _l_a));
  }
};

template <typename FA,typename FB>
class FunctionProduct 
{
private:
  FA _fn_a;
  FB _fn_b;

public:

  /**
   * copy constructor for function product
   **/
  void operator = (const FunctionProduct &f){
    _fn_a = f._fn_a;
    _fn_b = f._fn_b;
  };

  FunctionProduct(){};

  /**
   * constructor for product of functions
   **/
  FunctionProduct(FA fn_a, FB fn_b)
  {
    _fn_a = fn_a;
    _fn_b = fn_b;
  }

  /**
   * computes product of 2 functions fn_a and fn_b at x
   **/
  template <typename T>
  T operator()(T x)
  {

    T res =  _fn_a(x) * _fn_b(x);
    return res;
  }

  /**
   * first derivitive using the product rule
   **/
  template <typename T>
  T first_deriv(T x)
  {
    return _fn_a.first_deriv(x) * _fn_b(x) + _fn_b.first_deriv(x) * _fn_a(x);
  }

  /**
   * second derivative using the product rule: d^2/dx^2(f(x) g(x)) = g(x) f''(x) + 2 f'(x) g'(x) + f(x) g''(x)
   **/
  template <typename T>
  T second_deriv(T x)
  {
    return _fn_b(x) * _fn_a.second_deriv(x) + 2 * _fn_a.first_deriv(x) * _fn_b.first_deriv(x) + _fn_a(x) * _fn_b.second_deriv(x);
  }
};

template <typename F>
class Integral
{
private:
  F _fn;

public:

  Integral(){

  } 
  
  /**
   * constructor for integral
   **/
  Integral(F fn)
  {
    _fn = fn;
  }

  /*
   * calculates the integeral using the extended trapazoidal rule
   **/
  double operator()(int a, int b, int n)
  {
    double h = (1.0 * (b - a)) / (1.0 * n);
    double res = _fn(1.0 * a);
    for (int i = 1; i < n; i++)
    {
      res += (2.0 *_fn(a + i * h));
    }
    res += _fn(a + n * h);
    res /= (2.0 * n);
    res *= (b - a);
    return res;
  }

  /**
   * calcualtes truncation error at a given x based on the provided parameters
   **/
  double truncation_error(int a, int b, int n, double x)
  {
    return pow(b - a, 3) * _fn.second_deriv(x) / pow(n, 2);
  }
};


class Problem1TestCase{
  private:
    Integral <FunctionProduct<Gaussian,Gaussian>> _test_case;
    double _expected_result;
    int _id;

  public:
    Problem1TestCase(int id, double l_a, double l_b, double x_a, double x_b, double alpha, double beta, double expected_result){
      _id = id;
      _expected_result = expected_result;
      Gaussian g_a(x_a, l_a, alpha);
      Gaussian g_b(x_b, l_b, beta);
      FunctionProduct<Gaussian,Gaussian> fn(g_a, g_b);
      Integral<FunctionProduct<Gaussian,Gaussian>> test_case(fn);
      _test_case = test_case;
    }

    void operator()(int a, int b, int n){
      double result = _test_case(a, b, n); 
      std::cout << "\033[0;44m Test" << _id << " \33[0m " << std::setw(16) << result << " status: ";
      if(abs(result - _expected_result) < 0.0000001f){
	std::cout << "\033[32m PASSED" << "\033[0m" << std::endl;
      }
      else{
	std::cout << "\033[31m FAILED" << "\033[0m" << std::endl;
      }
    }
};

int main()
{

  // Problem One Tests
  std::cout << "Problem 1: " << std::endl;
  std::cout << "using a = -5 b = 5 and n = 10000 " <<std::endl;
  Problem1TestCase case1(1,0.0,0.0,0.0,0.0,1.0,1.0,1.25331413731550012);
  Problem1TestCase case2(2,0.0,1.0,0.0,0.0,1.0,1.0,0.00000000000000000);
  Problem1TestCase case3(3,0.0,0.0,0.0,1.0,1.0,1.0,7.60173450533140338e-01);
  Problem1TestCase case4(4,0.0,1.0,0.0,1.0,1.0,1.0,-3.80086725266570169e-01);
  case1(-5,5,10000);
  case2(-5,5,10000);
  case3(-5,5,10000);
  case4(-5,5,10000);

  // Problem 2 Tests
  std::cout << "Problem 2: " << std::endl;


}
