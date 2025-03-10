#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// SLC, Nov 9, 2011

#ifndef _CGOPTIMIZE_H_
#define _CGOPTIMIZE_H_

#define CGOPT_LINE_SEARCH_SECANT 0
#define CGOPT_LINE_SEARCH_BISECTION 1

#include "parstream.H"
#include "REAL.H"
#include <cmath>
#include "NamespaceHeader.H"
#include "CH_assert.H"

/** @file */ 
inline Real sqr(const Real& a)
{
  return a*a;
}


/// Specify operations required to carry out gradient based optimization
template <class X>
class ObjectiveWithGradient {

public:

  ///virtual destructor
   virtual ~ObjectiveWithGradient(){};

  /** 
      Compute dot product a.b
      @return a.b
  */
  virtual Real dotProduct(X& a , const X& b) = 0; 


  /**
     Compute a main objective (fm), a penalty function (fp) and the 
     gradient (g) of f + p with respect to x

     @param fm (output) main objective
     @param fp (output) penalty function
     @param g  (output) gradient
     @param x  (input)  state vector
     @param inner  (input)  is this an 'inner' iteration?
  */
   virtual void computeObjectiveAndGradient(Real& fm, Real& fp ,X& g, const X& x , bool inner) = 0;

  /// Copy x to y
   virtual void assign(X& y, const X& x) = 0;

  /// Set y = y + s * x
   virtual void incr(X& y, const X& x, Real s) = 0;

  /// resize a / allocate storage to conform with b
   virtual void create (X &a, const X &b) = 0;

  /// free storage 
   virtual void free (X &a) = 0;

  /// set a = a * s
   virtual void scale(X &a, const Real s) = 0;

  /// apply precoditioner K, so a = Kb
   virtual void preCond(X &a, const X &b) = 0;

  /// carry out whatver operations are required prior to a restart
   virtual void restart() = 0;

  /// number of degrees-of-freedom (maximum number of iterations before restart)
   virtual int nDoF(const X& x) = 0;

};

/// bisection method: find y = x + s*d such that |f'(y).f'(x)| << f'(x).f'(x)
/**
   x - initial point
   r - input - -f'(x) (steepest descent) 
   d - input - (preconditioned) search direction / steepest descent ( d = M^-1 r)
   y - output
   w - workspace
   a_fx - f(x) (where f' is the gradient of scalar f > 0)

 */
template < class X>
void bisectionLineSearch
(ObjectiveWithGradient<X>& a_F, X& x, X& r, X& d, X& y, X& w, Real a_fx, 
 int a_maxIter, Real a_initialStep, Real a_stepMaxGrow, Real a_tol)
{

  Real a, b, c, s; // brackets, mid-point, increment along d
  Real fa, fb, fc; // f'(x + s d).d @ s = a,b,c
  Real Fm, Fp; // misfit and penalty norms (F = Fm + Fp)
  Real f0, s0; // initial values

  int iter(0);
  a = 0.0;
  fa = a_F.dotProduct(r,d);
  f0 = fa;
  s0 = a_initialStep;
  s = s0;
  b = a;
  //walk the line until fb < 0 (since fa > 0) 
  do {
    b += s;
    a_F.assign(y, x);
    a_F.incr(y, d, b);
    a_F.computeObjectiveAndGradient(Fm,Fp,w,y,true);
    fb = -a_F.dotProduct(w,d);

    pout() << " ... bisection: searching for a bracket: iter " << iter
	   << " f'(x + bd).d = " << fb
	   << " f'(x).d = " << f0
	   << " b = " << b << " f(x + bd) = " << Fm + Fp << std::endl;
    
    if (fb < fa)
      {
	if (fb > 0.0)
	  {
	    // might as well move the lower limit 
	    fa = fb; a = b;	
	  }
      }
    else
      {
	// if d really is a descent direction, a reduction should help
        // need to start again?
	b = 0.0;
	a = 0.0;;
	s *= 0.5;
      }  
    
    iter++;
  } while ( (2*iter < a_maxIter) && (fb > 0.0) && (s > a_tol * s0 ) );
  
  do {
    c = 0.5*(a + b); //midpoint
    a_F.assign(y, x);
    a_F.incr(y, d, c);
    a_F.computeObjectiveAndGradient(Fm,Fp,w,y,true);

    fc = -a_F.dotProduct(w,d);
    
    pout() << " ... bisection iter " << iter 
	   << "  f'(x + sd).d = " << fc
	   << " f'(x).d = " << f0
	   << " s = " << c << " f(x + sd) = " << Fm + Fp << std::endl;
    
      if (std::signbit(fc) == std::signbit(fa))
      {
	a = c;
	fa = fc;
      }
    else
      {
	b = c;
	fb = fc;
      }
    iter++;
    
  } while ( (iter < a_maxIter) && (Abs(fc) > sqr(a_tol) * Abs(f0)) && (Abs(b-a) > a_tol * s0));
  a_F.assign(x, y);

}



/// Secant method: find y = x + s*d such that |f'(y).f'(x)| << f'(x).f'(x)
/**
   x - initial point
   r - input - ||f'(x)||^2
   d - input - search direction
   y - output
   p - workspace
   a_fx - f(x) (where f' is the gradient of scalar f > 0)

 */
template < class X>
void secantLineSearch
(ObjectiveWithGradient<X>& a_F, X& x, X& r, X& d, X& y, X& p, Real a_fx, 
 int a_maxIter, Real a_initialStep, Real a_stepMaxGrow, Real a_tol)
{

  //secant line search
  Real deltaD = a_F.dotProduct(r,d);
  Real eta = -deltaD;
  Real initialStep = a_initialStep;
  Real alpha = -initialStep;
  a_F.assign(y,x);
  a_F.incr(y,d,initialStep);
  Real fms,fps;
  a_F.computeObjectiveAndGradient(fms,fps,p,y,true);
  Real etaPrev = a_F.dotProduct(p,d);
  
  pout() << " ... initial secant ||f'(x+sd).f'(x)||^2 = " 
	 << -etaPrev << " s = " << initialStep << " f(x+sd) = " << fms << std::endl;
  
  while ( (fms + fps <  a_fx) && (-etaPrev > -eta))
    {
      //we have a problem, because although we have a descent direction in the  objective
      //function f , we have an ascent direction for f'(x+sd).f'(x). For now, just walk along d till
      // we do have descent in both, then start the secant method offset by alpha + a_initialStep
      pout() << " d is a descent direction for f(x+sd) but not for f'(x+sd).f'(x)  " 
	     << " so walking the line (remind SLC to look up a better method) " << std::endl;
      alpha -= initialStep;
      a_F.incr(y,d,initialStep);
      a_F.computeObjectiveAndGradient(fms,fps,p,y,true);
      eta = etaPrev;
      etaPrev = a_F.dotProduct(p,d);
      pout() << " ... initial secant f'(x+ sd).f'(x) = " 
	     << -etaPrev << " s = " << alpha << " f(x+ sd) = " << a_fx << std::endl;
    }
  
  while ( (fms + fps >  a_fx) && initialStep > 0.015624 * a_initialStep  )
    {
      // here we have an ascent direction the objective function. Try reducing the step.
      // this can of course be a happy failure
      pout() << " d is a not descent direction for f(x+sd), reducing s ";
      a_F.incr(y,d,-0.5*initialStep);
      initialStep *= 0.5;
      alpha = - initialStep;
      a_F.computeObjectiveAndGradient(fms,fps,p,y,true);
      eta = etaPrev;
      etaPrev = a_F.dotProduct(p,d);
      pout() << " ... initial secant f'(x+ sd).f'(x) = " 
	     << -etaPrev << " s = " << initialStep << " f(x+ sd) = " << a_fx << std::endl;
    }
  
  int j = 0;
  Real q = -(alpha + a_initialStep); // normally alpha = - a_initialStep, but sometimes we needed to walk the line
  if (q > 0.0)
    a_F.incr(x,d,q);
  
  alpha = -a_initialStep;
  
  do {
    
    if (j > 0)
      {
	a_F.computeObjectiveAndGradient(fms,fps,r,x,true);
	eta = a_F.dotProduct(r,d);
      }
    
    pout() <<  " ... secant iteration j = " << j << " f'(x[0]+qd).f'(x[0]) = " << -eta
	   << " q = " << q 
	   << std::endl;
	
    Real absAlphaPrev = std::abs(alpha);
    alpha = alpha * eta / (etaPrev-eta);
    
    //limit the rate at which alpha grows
    if ( (a_stepMaxGrow > 0.0) && (std::abs(alpha) > a_stepMaxGrow * absAlphaPrev) )
      {
	alpha = ((alpha > 0.0)?(1.0):(-1.0)) * a_stepMaxGrow * absAlphaPrev;
      }
    
    q += alpha;
    pout() << " ... secant step a = " << alpha << "  q = " << q  << std::endl;
    
    a_F.incr(x,d,alpha);
    etaPrev = eta;
    j++;
  } while (j < a_maxIter  && (sqr(alpha) * deltaD > sqr(a_tol)));
  
  
}



/**
   Nonlinear Conjugate-gradient method, used to mimimise F(x) w.r.t x 
*/
template <class X>
int CGOptimize(ObjectiveWithGradient<X>& a_F, X& a_x, 
	       int  a_maxIter,
	       Real a_tol,
	       Real a_hang,
	       int a_lineSearchMethod,
	       Real a_lineSearchInitialStep,
	       Real a_lineSearchStepMaxGrow,
	       int  a_lineSearchMaxIter,
	       Real a_lineSearchTol,
	       int a_iter = 0)
{

  //pout() << "CGOptimize -- a_iter = " << a_iter << endl;
  //Preconditioned Nonlinear Conjugate Gradients 
  //with Secant line OR bisection search and Polak-Ribiere updates
  //from J. R. Shewchuk, 1994 
  //"An Introduction to the Conjugate Gradient Method
  // Without the Agonizing Pain"

  X& x = a_x;
  X r,s,d,y,p;
  a_F.create(r,x);
  a_F.create(s,x);
  a_F.create(d,x);
  a_F.create(y,x);
  a_F.create(p,x);

  Real fm,fp, fmOld; // just for retreiving the the misfit norm and penalty norm
  a_F.computeObjectiveAndGradient(fm,fp,r,x,false);
  fmOld = fm;
  a_F.scale(r,-1.0);
  a_F.preCond(s,r);
  a_F.assign(d,s);
  Real deltaNew = a_F.dotProduct(r,s);
  Real deltaZero = deltaNew;
  int iter = a_iter; int k = 0;
  bool hang(false);
  // if we have a_iter > a_maxIter, we have presumably restarted
  int niter = a_maxIter;
  if ( (iter > a_maxIter) && (a_maxIter > 0))  niter += (iter/a_maxIter) * a_maxIter;
  
  while ( (iter < niter)
	  && (deltaNew > sqr(a_tol) * deltaZero) 
	  && (fm > 1.0e-16) 
	  && (!hang))
    {
      
      pout() << "CGOptimize iteration " << iter 
	     << " ||f'(x)||^2 = " << deltaNew 
	     << " ||fm(x)||^2 = " << fm 
	     << " ||fp(x)||^2 = " << fp 
	     << " ||fm(x)||^2 + ||fp(x)||^ = " << fm + fp;

      if (iter - a_iter > 0)
	{
	  pout() << " ||fm(x)||^2/||fm_old(x)||^2 = " << fm/fmOld;
	  hang = (fm > a_hang * fmOld);
	}
      pout() << std::endl;

      fmOld = fm;
      //fpOld = fp;

      //find the (approximate) root of f'(y).f'(x)
      if (a_lineSearchMethod == CGOPT_LINE_SEARCH_BISECTION)
	{
	  bisectionLineSearch(a_F, x, r, d, y, p, fm+fp, a_lineSearchMaxIter, 
			      a_lineSearchInitialStep, a_lineSearchStepMaxGrow, a_lineSearchTol);
	}
      else
	{
	  // we had secant as the only option since 2011, so make it default
	  secantLineSearch(a_F, x, r, d, y, p, fm+fp,  a_lineSearchMaxIter, 
			   a_lineSearchInitialStep, a_lineSearchStepMaxGrow, a_lineSearchTol);
	}
	
      a_F.computeObjectiveAndGradient(fm,fp,r,x,false);
      a_F.scale(r,-1.0);
      Real deltaOld = deltaNew;
      Real deltaMid = a_F.dotProduct(r,s);
      a_F.preCond(s,r);
      deltaNew = a_F.dotProduct(r,s);
      Real beta = (deltaNew - deltaMid)/deltaOld;
      k++;
      if (k == a_F.nDoF(x) || beta <= 0.0 )
	{
	  pout() << "CGOptimize restart k = " << k 
		 << " beta = " << beta << std::endl;
	  
	  a_F.assign(d,s);
	  a_F.restart();
	  k = 0;
	}
      else
	{
	  a_F.scale(d,beta);
	  a_F.incr(d,s,1.0);
	}
      iter++;
    }
  
  pout() << "CGOptimize iteration " << iter 
	 << " ||f'(x)||^2 = " << deltaNew 
	 << " ||fm(x)||^2 = " << fm 
	 << " ||fp(x)||^2 = " << fp   
	 << " ||fm(x)||^2 + ||fp(x)||^ = " << fm + fp
	 << " ||fm(x)||^2/||fm_old(x)||^2 = " << fm/fmOld 
	 << std::endl;
  

  a_F.free(r);
  a_F.free(s);
  a_F.free(d);
  a_F.free(y);
  a_F.free(p);

  return iter;

}


#include "NamespaceFooter.H"
#endif /*_CGOPTIM_H_*/
