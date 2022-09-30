#import <armadillo>
#import <iostream>
#import <cmath>
#import <stdio.h>
using arma::vec;
using arma::mat;
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

  template <typename T> T first_deriv(T x){
    return _fn(x);
  }

  template <typename T> T second_deriv(T x){
    return _fn.second_deriv(x);
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

class Shell{
  private:
    vec _r_a; //center (x,y,z, ...)
    mat _l_a; //matrix containing all possible (l,m,n....) combinations
    double _alpha;
  public:
    Shell(){};

    Shell(vec r_a, mat l_a, double alpha){
      _r_a = r_a;
      _l_a = l_a;
      _alpha = alpha;
    }

    mat operator ()(vec r){
      mat term1(_r_a.n_elem, _l_a.n_rows);
      for(int i = 0; i < _l_a.n_rows; i++){
	term1(i) = prod(pow(r - _r_a, _l_a(i)));
      }
      vec term2 = exp((-1.0 * _alpha) * pow((r-_r_a),2));
      return term1 * term2;
    }

    vec r_a(){
      return _r_a;
    }

    mat l_a(){
      return _l_a;
    }

    double alpha(){
      return _alpha;
    }

};

//util methods 

// calcualtes the produc for all integers from one to n
int factorial(int n){
  int product = n;
  for(int i = n-1; i > 0; i--){
    product *= i;
  }
  if(product < 1){
    return 1;
  }
  return product;
}

//calcualtes the a factorial where the output is the product of all integers from 1 to n where x%a == n%a
int factorial(int n, int a){
  int product = n;
  for(int i = n; i > 0; i-=a){
    product *= i;
  }
  if(product < 1){
    return 1;
  }
  return product;
} 


struct bar_is_space : std::ctype<char> {
  bar_is_space() : std::ctype<char>(get_data()) {}
  static mask const* get_data()
  {
    static mask rc[table_size];
    rc['|'] = std::ctype_base::space;
    return &rc[0];
  }
};

class ShellOverlapIntegral{
  private:
    Shell _s_a;
    Shell _s_b;


    double alphas_product(){
      return _s_a.alpha() * _s_b.alpha();
    }

    double alphas_sum(){
      return _s_a.alpha() + _s_b.alpha();
    }

    double dim_dist_sqr(int dim){
      return pow(_s_a.r_a()(dim)-_s_b.r_a()(dim),2); 
    }

    double exponential_prefactor(int dim){
      return exp(-1.0 * ( (alphas_product() * dim_dist_sqr(dim)) / alphas_sum()));
    }

    double root_term(){
      return sqrt(M_PI/alphas_sum());
    }

    double calc_binomial(int m, int n){
      return 1.0*factorial(m)/(1.0*factorial(n)*factorial(m-n));
    }

    double overlap_summation(double x_p, int l_pair_a, int l_pair_b, int dim){
      double summation = 0.0;
      for(int i = 0; i <= l_pair_a; i++){
	for(int j = 0; j <= l_pair_b; j++){
	  if((i+j)%2 == 0){
	    double binomial_term = calc_binomial(l_pair_a,i)*calc_binomial(l_pair_b,j);
	    double factorial_term = factorial(i+j-1,2);
	    double a_term = pow(x_p - _s_a.r_a()(dim), l_pair_a-i);
	    double b_term = pow(x_p - _s_b.r_a()(dim), l_pair_b-j);
	    double denominator = pow(2.0*alphas_sum(),(i+j)/2.0);
	    double step = (binomial_term*((factorial_term*a_term*b_term)/denominator));
	    summation+=step;
	  }
	}
      }
      return summation;
    }

    double product_center(int dim){
      return (_s_a.r_a()(dim)*_s_a.alpha()+_s_b.r_a()(dim)*_s_b.alpha())/alphas_sum();
    }
  
  public:
    ShellOverlapIntegral(Shell s_a, Shell s_b){
      _s_a = s_a;
      _s_b = s_b;
    };

    mat operator ()(){
      mat result(_s_a.l_a().n_rows, _s_b.l_a().n_rows, arma::fill::ones);
      for(int i = 0; i < _s_a.r_a().n_elem; i++){
	double x_p = product_center(i);
	for(int k = 0; k < _s_a.l_a().n_rows; k++){
	  for(int l = 0; l < _s_b.l_a().n_rows; l++){ 
	    double z = exponential_prefactor(i) * root_term() * overlap_summation(x_p, _s_a.l_a()(k,i),_s_b.l_a()(l,i),i);
	    result(k,l) *= z;    
	  }
	}
      }
      return result;
    }
};
/**
 * @brief reads files with shell data
 *
 * @detail reads molecule files where the first line is the number of shells and
 *         the remaining lines follow the format below:
 *         x y z ...|l_a0 l_a1 l_a3 ...;...|alpha
 *
 *
 *
 *         the number of shells described by the file must match the specified number of shells in the first line
 * 
 **/
void read_shells(string &filename) {
  std::ifstream infile(filename);
  infile.imbue(std::locale(std::cin.getloc(), new bar_is_space));
  if (infile.is_open()) {
    std::string r_a_str;
    std::string l_a_str;
    double alpha;
    //read first line
    int expectedShells;
    infile >> expectedShells;;
    std::cout << "expecting " << expectedShells << " shells" << std::endl;
    int shellsCounted = 0;
    infile >> r_a_str >> l_a_str >> alpha;
    vec r_a(r_a_str);
    mat l_a(l_a_str);
    Shell s_a(r_a,l_a,alpha);
    std::string r_b_str;
    std::string l_b_str;
    double beta;
    infile >> r_b_str >> l_b_str >> beta;
    vec r_b(r_b_str);
    mat l_b(l_b_str);
    Shell s_b(r_b,l_b,beta);

    ShellOverlapIntegral overlap(s_a,s_b);
    std::cout << overlap() << std::endl;
    infile.close();
  } 
  else {
    throw std::invalid_argument("Can't open file to read.");
  }
}

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
  std::string filename = "shelltest.in";
  read_shells(filename);
  return 0;
}
