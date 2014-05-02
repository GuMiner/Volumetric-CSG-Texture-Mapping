#pragma once
/**************************************************************************/
/* File:   opti.hpp                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

namespace netgen
{

  /** 
      Function to be minimized.
  */
  class MinFunction 
  {
  public:
    ///
    virtual double Func (const Vector & x) const;
    ///
    virtual void Grad (const Vector & x, Vector & g) const;
    /// function and gradient
    virtual double FuncGrad (const Vector & x, Vector & g) const;
    /// directional derivative
    virtual double FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const;
    /// if |g| < gradaccuray, then stop bfgs
    virtual double GradStopping (const Vector & /* x */) const { return 0; } 

    ///
    virtual void ApproximateHesse (const Vector & /* x */,
				   DenseMatrix & /* hesse */) const;
  };

  class OptiParameters
  {
  public:
    int maxit_linsearch;
    int maxit_bfgs;
    double typf;
    double typx;

    OptiParameters ()
    {
      maxit_linsearch = 100;
      maxit_bfgs = 100;
      typf = 1;
      typx = 1;
    }
  };
  
  /** Implementation of BFGS method.
      Efficient method for non-linear minimiztion problems.
      @param x initial value and solution 
      @param fun function to be minimized
  */
  extern double BFGS (Vector & x, const MinFunction & fun, 
		      const OptiParameters & par,
		      double eps = 1e-8);

  /** Steepest descent method.
      Simple method for non-linear minimization problems.
      @param x initial value and solution 
      @param fun function to be minimized
  */
  void SteepestDescent (Vector & x, const MinFunction & fun,
			const OptiParameters &  par);

  extern void lines (
		     Vector & x,         // i: Ausgangspunkt der Liniensuche
		     Vector & xneu,      // o: Loesung der Liniensuche bei Erfolg
		     Vector & p,         // i: Suchrichtung
		     double & f,         // i: Funktionswert an der Stelle x
		     // o: Funktionswert an der Stelle xneu, falls ifail = 0
		     Vector & g,         // i: Gradient an der Stelle x
		     // o: Gradient an der Stelle xneu, falls ifail = 0

		     const MinFunction & fun,  // function to minmize
		     const OptiParameters & par, // parameters
		     double & alphahat,  // i: Startwert f�r alpha_hat
		     // o: Loesung falls ifail = 0
		     double fmin,        // i: untere Schranke f�r f
		     double mu1,         // i: Parameter mu_1 aus Alg.2.1
		     double sigma,       // i: Parameter sigma aus Alg.2.1
		     double xi1,         // i: Parameter xi_1 aus Alg.2.1
		     double xi2,         // i: Parameter xi_1 aus Alg.2.1
		     double tau,         // i: Parameter tau aus Alg.2.1
		     double tau1,        // i: Parameter tau_1 aus Alg.2.1
		     double tau2,        // i: Parameter tau_2 aus Alg.2.1
		     int & ifail);        // o:  0 bei erfolgreicher Liniensuche
  //    -1 bei Abbruch wegen Unterschreiten von fmin
  //    1 bei Abbruch, aus sonstigen Gr�nden
}