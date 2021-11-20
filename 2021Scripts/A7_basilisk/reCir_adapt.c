/**
# Shock reflection by a circular cylinder

The evolution of an initial "step" wave is modelled using the
Saint-Venant equations. The wave interacts with a circular cylinder
described using embedded solid boundaries. Adaptivity is used to track
the wave fronts. This example is discussed in [An and Yu,
2012](/src/references.bib#an2012). */

#include "saint-venant.h"
/** gfs and vtk output */
#include "vtk.h"
/** don't know how this works. The program complains. */
//#include "output_vtk_foreach_new.h"

int LEVEL = 8;
int MAXLEVEL = 9;
int MINLEVEL = 5;

/**
We define a new boundary for the cylinder. */

// bid cylinder;

int main() {
  size (7.05);
  G = 9.81;
  // origin (-L0/2., -L0/2.);
  init_grid (1 << LEVEL);
  CFL = 0.490;
  theta = 1.5; // use Sweby limiter
  run();
}

/**
We impose height and velocity on the left boundary. */

#define H0 0.02 // 2cm depth
#define U0 0.20
#define INLETY 0.705
#define LY 1.41
#define CF 0.007

double simTime = 1000;

// h[left]   = y>INLETY?H0:;
// eta[left] = H0;
u.n[left] = y>INLETY?U0:0.0;
u.t[left] = 0.0;

u.n[right]  = + radiation(0);
// h[right] = neumann(0.);

scalar conc[];

conc[right] = neumann(0);

event init (i = 0) {

  /**
  The geometry is defined by masking and the initial step function is
  imposed. */
  mask(y > LY ? top : none);
  // mask (sq(x + 0.5) + sq(y) < sq(0.5) ? cylinder : none);
  foreach() {
    h[] = H0;
    u.x[] = 0.0;
    u.y[] = 0.0;
    conc[] = 1.0;
  }
}

event friction(i++)
{
  // double uMed;
  
  foreach ()
  {
    // double a = h[] < dry ? HUGE : 1. + (cf / (2.)) * dt * norm(u) / h[];
      double a = h[] < dry ? HUGE : 1. + (CF / (2.)) * dt * norm(u) / h[];
      u.x[] = (u.x[] ) / a;
      u.y[] /= a;
      // uMed = u.x[] + dt *  (-(cf / 2.) * u.x[] * norm(u) / h[] + gravityCoeff*So);
      // uMed = (3. / 4.) * u.x[] + (1. / 4.) * uMed + (1. / 4.) * dt * (-(cf / 2.) * u.x[] * sqrt(sq(uMed)+sq(u.y[])) / h[] + gravityCoeff*So);
      // u.x[] = (1. / 3.) * u.x[] + (2. / 3.) * uMed + (2. / 3.) * dt * (-(cf / 2.) * u.x[] * sqrt(sq(uMed)+sq(u.y[])) / h[] + gravityCoeff*So);    

      // uMed = u.y[] + dt *  (-(cf / 2.) * u.y[] * norm(u) / h[]);
      // uMed = (3. / 4.) * u.y[] + (1. / 4.) * uMed + (1. / 4.) * dt * (-(cf / 2.) * u.y[] * sqrt(sq(uMed)+sq(u.x[])) / h[]);
      // u.y[] = (1. / 3.) * u.y[] + (2. / 3.) * uMed + (2. / 3.) * dt * (-(cf / 2.) * u.y[] * sqrt(sq(uMed)+sq(u.x[])) / h[]);  
  }
}

/**
We generate movies of depth and level of refinement. */

// event movie (t += 0.5; t <= simTime) {
//   output_ppm (h, min = 0.50001, max = 1.25*H0, map = cool_warm, linear = true,
// 	      n = 600, file = "depth.mp4");
//   scalar l[];
//   foreach()
//     l[] = level;
//   output_ppm (l, map = cool_warm, min = 4, max = LEVEL, n = 400,
// 	      file = "level.mp4");
// }

//adaptivity
event adapt(i++)
{
  astats s = adapt_wavelet({h, u.x, u.y}, (double[]){H0 / 100.0, U0/100.0, U0/100.0}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}

// event adapt (i++) {
//  astats s = adapt_wavelet ({h}, (double[]){1e-2}, LEVEL);
//  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
//}

event outputVTKFiles (t += 100.00; t <= simTime) {
    char name1[80];
    char name2[80];
     sprintf(name1, "out-%.1f.txt", t);
    sprintf(name2, "out-%.1f.gfs", t);
    FILE *fp1 = fopen(name1, "w");
  foreach ()
  {
      fprintf(fp1, "%g %g %g %g %g\n", x, y, h[], u.x[], u.y[]);
  }
    FILE *fp2 = fopen(name2, "w");
    output_gfs(fp2, translate = true);
}

/**
## See also

* [Same case with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/shock.html)
*/
