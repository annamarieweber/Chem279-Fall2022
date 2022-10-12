#if !defined NUMERICAL_H
#define NUMERICAL_H

#define PI 3.1415926

class Integrand
{
public:
virtual double eval(double x) = 0;
};

double DEInfinity(Integrand& f1, Integrand& f2,const int nmupoints_oneside, const double stepsize);

class Gussian: public Integrand
{
    private:
    double x0;
    double alpha;
    int l;
    public:
    Gussian(double x0_input, double alpha_input, int l_input):
    x0(x0_input), alpha(alpha_input), l(l_input) {}
    Gussian():x0(0.0), alpha(0.5), l(0){}
    ~Gussian(){}
    double eval(double x);
};

#endif // NUMERICAL_H