#include "evaluation.h"
#include "cluster.h"
#include <sciplot/sciplot.hpp>
#include <armadillo>
#include <valarray>
using namespace sciplot;
using arma::mat;
using arma::vec;

double calcTruncationError(Cluster c, double h, int method)
{
  double e;
  if (method == 1)
  {
    e = norm(c.calcAnalyticalForce() - c.calcCentralFdForce(h));
  }
  else
  {
    e = norm(c.calcAnalyticalForce() - c.calcForwardFdForce(h));
  }
  std::cout << "e val" << std::endl;
  std::cout << e << std::endl;
  std::cout << log(e) << std::endl;
  return log(e);
}

void evaluateFDApproximation(Cluster c)
{
  struct truncationError
  {
    truncationError(){};
    truncationError(Cluster cluster, int m)
    {
      c = cluster;
      method = m;
    }

    std::vector<double> operator()(vec x) const
    {
      std::vector<double> output;
      for (int i = 0; i < x.n_elem; i++)
      {
        output.push_back(calcTruncationError(c, x(i), method));
      }
      return output;
    }

  private:
    Cluster c;
    int axis;
    int method;
  };

  std::vector<double> x = {0.1, 0.01, 0.001, 0.0001};
  std::vector<double> log_x = {log(0.1), log(0.01), log(0.001), log(0.0001)};

  // Create a Plot object
  Plot2D plot;

  // Set the x and y labels
  plot.ylabel("Log(TruncationError)");
  plot.xlabel("Log(stepsize)");

  // Set the x and y ranges
  plot.xrange(-5.0, -2.0);
  plot.yrange(0.0, 20.0);

  // Set the legend to be on the bottom along the horizontal
  plot.legend()
      .atOutsideBottom()
      .displayHorizontal()
      .displayExpandWidthBy(2);

  // Plot sin(i*x) from i = 1 to i = 6

  plot.drawCurve(log_x, truncationError(c, 0)(x)).label("Truncation error - ForwardFD");
  plot.drawCurve(log_x, truncationError(c, 1)(x)).label("Truncation error - CentralFD");

  // Create figure to hold plot
  Figure fig = {{plot}};
  // Create canvas to hold figure
  Canvas canvas = {{fig}};

  // Show the plot in a pop-up window
  canvas.show();

  // Save the plot to a PDF file
  canvas.save("../../output/evaluateFDApproximation.pdf");
}
