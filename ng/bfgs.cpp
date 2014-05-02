/***************************************************************************/
/*                                                                         */
/* Vorlesung Optimierung I, Gfrerer, WS94/95                               */
/* BFGS-Verfahren zur Lösung freier nichtlinearer Optimierungsprobleme     */
/*                                                                         */
/* Programmautor:  Joachim Schöberl                                        */
/* Matrikelnummer: 9155284                                                 */
/*                                                                         */
/***************************************************************************/

#include "stdafx.h" 
#include "mystdlib.h"
#include "myadt.hpp"

#include "linalg.hpp"
#include "opti.hpp"

namespace netgen
{

void MultLDLt (const DenseMatrix & l, const Vector & d, const Vector & g, Vector & p)
{
  double val;

  int n = l.Height();
  p = g;
  
  for (int i = 0; i < n; i++)
    {
      val = 0;
      for (int j = i; j < n; j++)
	val += p(j) * l(j, i);
      p(i) = val;
    }

  for (int i = 0; i < n; i++)
    p(i) *= d(i);

  for (int i = n-1; i >= 0; i--)
    {
      val = 0;
      for (int j = 0; j <= i; j++)
	val += p(j) * l(i, j);
      p(i) = val;
    }
}

void SolveLDLt (const DenseMatrix & l, const Vector & d, const Vector & g, Vector & p)
{
  double val;

  int n = l.Height();
  p = g;

  for (int i = 0; i < n; i++)
    {
      val = 0;
      for (int j = 0; j < i; j++)
	val += p(j) * l(i,j);
      p(i) -= val;
    }

  for (int i = 0; i < n; i++)
    p(i) /= d(i);
  
  for (int i = n-1; i >= 0; i--)
    {
      val = 0;
      for (int j = i+1; j < n; j++)
	val += p(j) * l(j, i);
      p(i) -= val;
    }
}

int LDLtUpdate (DenseMatrix & l, Vector & d, double a, const Vector & u)
{
  // Bemerkung: Es wird a aus R erlaubt
  // Rueckgabewert: 0 .. D bleibt positiv definit
  //                1 .. sonst

  int n = l.Height();

  Vector v(n);
  double t, told, xi;

  told = 1;
  v = u;

  for (int j = 1; j <= n; j++)
    {
      t = told + a * sqr (v(j-1)) / d(j-1);

      if (t <= 0) 
	{
	  cout << "update err, t = " << t << endl;
	  return 1;
	}

      xi = a * v(j-1) / (d(j-1) * t);

      d(j-1) *= t / told;

      for (int i = j + 1; i <= n; i++)
	{
	  v(i-1) -= v(j-1) * l.Elem(i, j);
	  l.Elem(i, j) += xi * v(i-1);
	}

      told = t;
    }

  return 0;
}

double BFGS (
	     Vector & x,         // i: Startwert
	     // o: Loesung, falls IFAIL = 0
	     const MinFunction & fun,
	     const OptiParameters & par,
	     double eps
	     )


{
  int n = x.Size();
  long it;
  char a1crit, a3acrit;


  Vector d(n), g(n), p(n), temp(n), bs(n), xneu(n), y(n), s(n), x0(n);
  DenseMatrix l(n);
  DenseMatrix hesse(n);

  double /* normg, */ alphahat, hd, fold;
  double a1, a2;
  const double mu1 = 0.1, sigma = 0.1, xi1 = 1, xi2 = 10;
  const double tau = 0.1, tau1 = 0.1, tau2 = 0.6;

  Vector typx(x.Size());      // i: typische Groessenordnung der Komponenten
  double f, f0;
  double typf;               // i: typische Groessenordnung der Loesung
  double fmin = -1e5;           // i: untere Schranke fuer Funktionswert
  //  double eps = 1e-8;            // i: Abbruchschranke fuer relativen Gradienten
  double tauf = 0.1;            // i: Abbruchschranke fuer die relative Aenderung der
                                //    Funktionswerte
  int ifail;                    // o:  0 .. Erfolg
                                //    -1 .. Unterschreitung von fmin
                                //     1 .. kein Erfolg bei Liniensuche
                                //     2 .. Überschreitung von itmax

  typx = par.typx;
  typf = par.typf;


  l = 0;
  for (int i = 1; i <= n; i++)
    l.Elem(i, i) = 1;

  f = fun.FuncGrad (x, g);
  f0 = f;
  x0 = x;

  it = 0;
  do
    {
      // Restart
      if (it % (5 * n) == 0)
	{

	  for (int i = 1; i <= n; i++)
	    d(i-1) = typf/ sqr (typx(i-1));   // 1;
	  for (int i = 2; i <= n; i++)
	    for (int j = 1; j < i; j++)
	      l.Elem(i, j) = 0;

	}

      it++;
      if (it > par.maxit_bfgs)
	{
	  ifail = 2;
	  break;
	}


      // Solve with factorized B
      SolveLDLt (l, d, g, p);

      p *= -1;
      y = g;

      fold = f;

      // line search
      alphahat = 1;
      lines (x, xneu, p, f, g, fun, par, alphahat, fmin,
	     mu1, sigma, xi1, xi2, tau, tau1, tau2, ifail);

      if(ifail == 1)
	cout << "no success with linesearch" << endl;

      s.Set2 (1, xneu, -1, x);
      y *= -1;
      y.Add (1,g); // y += g;

      x = xneu;

      // BFGS Update
      MultLDLt (l, d, s, bs);

      a1 = y * s;
      a2 = s * bs;

      if (a1 > 0 && a2 > 0)
	{
	  if (LDLtUpdate (l, d, 1 / a1, y) != 0)
	    {
              cerr << "BFGS update error1" << endl;
	      cout << "BFGS update error1" << endl;
	      cout << "l " << endl << l << endl
			 << "d " << d << endl;
	      ifail = 1;
	      break;
	    }

	  if (LDLtUpdate (l, d, -1 / a2, bs) != 0)
	    {
              cerr << "BFGS update error2" << endl;
	      cout << "BFGS update error2" << endl;
	      cout << "l " << endl << l << endl
			 << "d " << d << endl;
	      ifail = 1;
	      break;
	    }
	}

      // Calculate stop conditions
      hd = eps * max2 (typf, fabs (f));
      a1crit = 1;
      for (int i = 1; i <= n; i++)
	if ( fabs (g(i-1)) * max2 (typx(i-1), fabs (x(i-1))) > hd)
	  a1crit = 0;


      a3acrit = (fold - f <= tauf * max2 (typf, fabs (f)));

      if (g.L2Norm() < fun.GradStopping (x)) break;

    }
  while (!a1crit || !a3acrit);


  if (f0 < f || (ifail == 1))
    {
      cout << "fail, f = " << f << " f0 = " << f0 << endl;
      f = f0;
      x = x0;
    }

  return f;
}
}
