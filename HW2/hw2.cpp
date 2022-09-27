#import <armadillo>
#import <stdio>
#import <cmath>

class Function{
  template <typename T> T operator () (T x){
    return 1 * x;
  }

  template <typename T> T first_deriv(T x){
    return 1 * x;
  }

  template <typename T> T second_deriv(T x){
    return 1 * x;
  }
};

class Gaussian: public Function {
  private:
    double _x_a;
    double _l_a;
    double _alpha;
  public:
    /**
     * constructor for a gaussian
     **/
    Gaussian(double x_a, double l_a, double alpha){
      _x_a = x_a;
      _l_a = l_a;
      _alpha = alpha;
    }

    /**
     * evaluation of the Gaussian Function = (x − XA)
     **/
    template <typename T> T operator () (T x){
      return pow(x- _x_a, _l_a) * exp(-1.0 * _alpha * pow(x - _x_a,2));
    };

    /**
     * first derivative of the gaussian
     * e^(-alpha (x - _x_a)^2) (x - _x_a)^(_l_a - 1) (_l_a - 2 alpha (x - _x_a)^2)
     **/
    template <typename T> T first_deriv(T x){
      return exp(-1.0 * _alpha * pow(x - _x_a, 2)) * pow(x - _x_a, _l_a - 1) * (_l_a - 2 * _alpha * pow(x - _x_a, 2));
    }

    /**
     * second derivative of the gaussian
     * e^(-_alpha (-_x_a + x)^2) (-_x_a + x)^(_l_a - 2) (4 _alpha^2 (_x_a - x)^4 - 2 _alpha (2 _l_a + 1) (-_x_a + x)^2 + _l_a (-1 + _l_a))
     **/
    template <typename T> T second_deriv(T x){
      return exp(-1.0 *_alpha * pow(-1.0 *(_x_a + x), 2) * pow(-1.0 * _x_a + x, _l_a - 2)) * (4 * pow(_alpha, 2) * pow(_x_a - x, 4) - 2 * _alpha * (2 * _l_a + 1)  * (-1 * pow(_x_a + x, 2)) + _l_a * (-1 + _l_a));
    }
};

class FunctionProduct: public Function{
  private:
    Function _fn_a;
    Function _fn_b;

  public:
    /**
     * constructor for product of functions
     **/
    FunctionProduct(Function fn_a, Function fn_b){
      _fn_a = fn_a;
      _fn_b = fn_b;
    }

    /**
     * computes product of 2 functions fn_a and fn_b at x
     **/
    template <typename T> T operator () (T x){
      return _fn_a(x) * _fn_b(x);
    }

    /**
     * first derivitive using the product rule
     **/
    template <typename T> T first_deriv(T x){
      return _fn_a.first_deriv(x) * _fn_b(x) + _fn_b.first_deriv(x) * _fn_a(x);
    }

    /**
     * second derivative using the product rule: d^2/dx^2(f(x) g(x)) = g(x) f''(x) + 2 f'(x) g'(x) + f(x) g''(x)
     **/
    template <typename T> T second_deriv(T x){
      return _fn_b(x) * _fn_a.second_deriv(x) + 2 * _fn_a.first_deriv(x) * _fn_b.first_deriv(x) + _fn_a(x) * _fn_b.second_deriv(x);
    }
};

int main(){

}
