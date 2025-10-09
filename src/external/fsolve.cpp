# include <cmath>
# include <cstdlib>
# include <iomanip>
# include <iostream>
# include <limits>

using namespace std;

# include "fsolve.hpp"

//****************************************************************************80

void backward_euler_residual ( void dydt ( double t, double y[], double dy[] ), 
  int n, double to, double yo[], double tm, double ym[], double fm[] )

//****************************************************************************80
//
//  Purpose:
//
//    backward_euler_residual() evaluates the backward Euler residual.
//
//  Discussion:
//
//    Let to and tm be two times, with yo and ym the associated ODE
//    solution values there.  If ym satisfies the backward Euler condition,
//    then
//
//      dydt(tm,ym) = ( ym - yo ) / ( tm - to )
//
//    This can be rewritten as
//
//      residual = ym - yo - ( tm - to ) * dydt(tm,ym)
//
//    Given the other information, a nonlinear equation solver can be used
//    to estimate the value ym that makes the residual zero.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    08 November 2023
//
//  Author:
//
//    John Burkardt
//
//  Input:
//
//    void dydt( double t, double y[], double dy[] ): 
//    evaluates the right hand side of the ODE.
//
//    int n: the vector size.
//
//    double to, yo[n]: the old time and solution.
//
//    double tm, ym[n]: the new time and tentative solution.
//
//  Output:
//
//    double fm[n]: the backward Euler residual.
//
{
  int i;
  double *dydtm;
  
  dydtm = new double[n];

  dydt ( tm, ym, dydtm );

  for ( i = 0; i < n; i++ )
  {
    fm[i] = ym[i] - yo[i] - ( tm - to ) * dydtm[i];
  }

  delete [] dydtm;

  return;
}
//****************************************************************************80

void bdf2_residual ( void dydt ( double t, double y[], double dy[] ), 
  int n, double t1, double y1[], double t2, double y2[], double t3, 
  double y3[], double fm[] )

//****************************************************************************80
//
//  Purpose:
//
//    bdf2_residual() evaluates the BDF2 residual.
//
//  Discussion:
//
//    Let t1, t2 and t3 be three times, with y1, y2 and y3 the associated ODE
//    solution values there.  Assume only the y3 value may be varied.
//
//    The BDF2 condition is:
//
//      w = ( t3 - t2 ) / ( t2 - t1 )
//      b = ( 1 + w )^2 / ( 1 + 2 w )
//      c = w^2 / ( 1 + 2 w )
//      d = ( 1 + w ) / ( 1 + 2 w )
//
//      y3 - b y2 + c y1 = d * ( t3 - t2 ) * dydt( t3, y3 )
//
//    but if (t3-t2) = (t2-t1), we have:
//
//      w = 1
//      b = 4/3
//      c = 1/3
//      d = 2/3
//      y3 - 4/3 y2 + 1/3 y1 = 2 dt * dydt( t3, y3 )
//
//    This can be rewritten as
//
//      residual = y3 - b y2 + c y1 - d * ( t3 - t2 ) * dydt(t3,y3)
//
//    This is the BDF2 residual to be evaluated.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    18 November 2023
//
//  Author:
//
//    John Burkardt
//
//  Input:
//
//    void dydt( double t, double y[], double dy[] ): 
//    evaluates the right hand side of the ODE.
//
//    int n: the vector size.
//
//    double t1, y1[n], t2, y2[n], t3, y3[n]: three sets of
//    data at a sequence of times.
//
//  Output:
//
//    double fm[n]: the backward Euler residual.
//
{
  double b;
  double c;
  double d;
  double *dydt3;
  int i;
  double w;

  w = ( t3 - t2 ) / ( t2 - t1 );
  b = ( 1.0 + w ) * ( 1.0 + w ) / ( 1.0 + 2.0 * w );
  c = w * w / ( 1.0 + 2.0 * w );
  d = ( 1.0 + w ) / ( 1.0 + 2.0 * w );

  dydt3 = new double[n];

  dydt ( t3, y3, dydt3 );

  for ( i = 0; i < n; i++ )
  {
    fm[i] = y3[i] - b * y2[i] + c * y1[i] - d * ( t3 - t2 ) * dydt3[i];
  }

  delete [] dydt3;

  return;
}
//****************************************************************************80

void dogleg ( int n, double r[], int lr, double diag[], double qtb[],
  double delta, double x[], double wa1[], double wa2[] )

//****************************************************************************80
//
//  Purpose:
//
//    dogleg() combines Gauss-Newton and gradient for a minimizing step.
//
//  Discussion:
//
//    Given an M by N matrix A, an n by n nonsingular diagonal
//    matrix d, an m-vector b, and a positive number delta, the
//    problem is to determine the convex combination x of the
//    gauss-newton and scaled gradient directions that minimizes
//    (a*x - b) in the least squares sense, subject to the
//    restriction that the euclidean norm of d*x be at most delta.
//
//    This function completes the solution of the problem
//    if it is provided with the necessary information from the
//    qr factorization of a. 
//
//    That is, if a = q*r, where q has orthogonal columns and r is an upper 
//    triangular matrix, then dogleg expects the full upper triangle of r and
//    the first n components of Q'*b.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    05 April 2010
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    This version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Parameters:
//
//    Input, int N, the order of R.
//
//    Input, double R[LR], the upper triangular matrix R stored by rows.
//
//    Input, int LR, the size of the storage for R, which should be at
//    least (n*(n+1))/2.
//
//    Input, double DIAG[N], the diagonal elements of the matrix D.
//
//    Input, double QTB[N], the first n elements of the vector 
//    (q transpose)*b.
//
//    Input, double DELTA, an upper bound on the euclidean norm of d*x.
//
//    Output, double X[N], contains the desired convex combination of the 
//    gauss-newton direction and the scaled gradient direction.
//
//    Workspace, WA1[N].
//
//    Workspace, WA2[N].
//
{
  double alpha;
  double bnorm;
  double epsmch;
  double gnorm;
  int i;
  int j;
  int jj;
  int jp1;
  int k;
  int l;
  double qnorm;
  double sgnorm;
  double sum;
  double temp;
//
//  EPSMCH is the machine precision.
//
  epsmch = std::numeric_limits<double>::epsilon();
//
//  Calculate the Gauss-Newton direction.
//
  jj = ( n * ( n + 1 ) ) / 2 + 1;

  for ( k = 1; k <= n; k++ )
  {
    j = n - k + 1;
    jp1 = j + 1;
    jj = jj - k;
    l = jj + 1;
    sum = 0.0;
    for ( i = jp1; i <= n; i++ )
    {
      sum = sum + r[l-1] * x[i-1];
      l = l + 1;
    }
    temp = r[jj-1];
    if ( temp == 0.0 )
    {
      l = j;
      for ( i = 1; i <= j; i++ )
      {
        temp = fmax ( temp, fabs ( r[l-1] ) );
        l = l + n - i;
      }
      temp = epsmch * temp;
      if ( temp == 0.0 )
      {
        temp = epsmch;
      }
    }
    x[j-1] = ( qtb[j-1] - sum ) / temp;
  }
//
//  Test whether the Gauss-Newton direction is acceptable.
//
  for ( j = 0; j < n; j++ )
  {
    wa1[j] = 0.0;
    wa2[j] = diag[j] * x[j];
  }
  qnorm = enorm ( n, wa2 );

  if ( qnorm <= delta )
  {
    return;
  }
//
//  The Gauss-Newton direction is not acceptable.
//  Calculate the scaled gradient direction.
//
  l = 0;
  for ( j = 0; j < n; j++ )
  {
    temp = qtb[j];
    for ( i = j; i < n; i++ )
    {
      wa1[i-1] = wa1[i-1] + r[l-1] * temp;
      l = l + 1;
    }
    wa1[j] = wa1[j] / diag[j];
  }
//
//  Calculate the norm of the scaled gradient and test for
//  the special case in which the scaled gradient is zero.
//
  gnorm = enorm ( n, wa1 );
  sgnorm = 0.0;
  alpha = delta / qnorm;
//
//  Calculate the point along the scaled gradient
//  at which the quadratic is minimized.
//
  if ( gnorm != 0.0 )
  {
    for ( j = 0; j < n; j++ )
    {
      wa1[j] = ( wa1[j] / gnorm ) / diag[j];
    }
    l = 0;
    for ( j = 0; j < n; j++ )
    {
      sum = 0.0;
      for ( i = j; i < n; i++ )
      {
        sum = sum + r[l] * wa1[i];
        l = l + 1;
      }
      wa2[j] = sum;
    }
    temp = enorm ( n, wa2 );
    sgnorm = ( gnorm / temp ) / temp;
    alpha = 0.0;
//
//  If the scaled gradient direction is not acceptable,
//  calculate the point along the dogleg at which the quadratic is minimized.
//
    if ( sgnorm < delta)
    {
      bnorm = enorm ( n, qtb );
      temp = ( bnorm / gnorm ) * ( bnorm / qnorm ) * ( sgnorm / delta );
      temp = temp - ( delta / qnorm ) * ( sgnorm / delta ) * ( sgnorm / delta )
        + sqrt ( pow ( temp - ( delta / qnorm ), 2 )
        + ( 1.0 - ( delta / qnorm ) * ( delta / qnorm ) )
        * ( 1.0 - ( sgnorm / delta ) * ( sgnorm / delta ) ) );
      alpha = ( ( delta / qnorm )
        * ( 1.0 - ( sgnorm / delta ) * ( sgnorm / delta ) ) ) / temp;
    }
  }
//
//  Form appropriate convex combination of the Gauss-Newton
//  direction and the scaled gradient direction.
//
  temp = ( 1.0 - alpha ) * fmin ( sgnorm, delta );
  for ( j = 0; j < n; j++ )
  {
    x[j] = temp * wa1[j] + alpha * x[j];
  }
  return;
}
//****************************************************************************80

double enorm ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    enorm() returns the Euclidean norm of a vector.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    05 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Input:
//
//    int N, the number of entries in A.
//
//    double X[N], the vector whose norm is desired.
//
//  Output:
//
//    double ENORM, the norm of X.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + x[i] * x[i];
  }
  value = sqrt ( value );
  return value;
}
//****************************************************************************80

void fdjac1 ( void fcn ( int n, double x[], double f[] ),
  int n, double x[], double fvec[], double fjac[], int ldfjac,
  int ml, int mu, double epsfcn, double wa1[], double wa2[] )

//****************************************************************************80
//
//  Purpose:
//
//    fdjac1() estimates an N by N Jacobian matrix using forward differences.
//
//  Discussion:
//
//    This function computes a forward-difference approximation
//    to the N by N jacobian matrix associated with a specified
//    problem of N functions in N variables. 
//
//    If the jacobian has a banded form, then function evaluations are saved 
//    by only approximating the nonzero terms.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    05 April 2010
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    This version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Parameters:
//
//    Input, void FCN(int n,double x[], double fx[]): the name of the 
//    C++ routine which returns in fx[] the function value at the n-dimensional
//    vector x[].
//
//    Input, int N, the number of functions and variables.
//
//    Input, double X[N], the evaluation point.
//
//    Input, double FVEC[N], the functions evaluated at X.
//
//    Output, double FJAC[N*N], the approximate jacobian matrix at X.
//
//       ldfjac is a positive integer input variable not less than n
//         which specifies the leading dimension of the array fjac.
//
//       ml is a nonnegative integer input variable which specifies
//         the number of subdiagonals within the band of the
//         jacobian matrix. if the jacobian is not banded, set
//         ml to at least n - 1.
//
//       epsfcn is an input variable used in determining a suitable
//         step length for the forward-difference approximation. this
//         approximation assumes that the relative errors in the
//         functions are of the order of epsfcn. if epsfcn is less
//         than the machine precision, it is assumed that the relative
//         errors in the functions are of the order of the machine
//         precision.
//
//       mu is a nonnegative integer input variable which specifies
//         the number of superdiagonals within the band of the
//         jacobian matrix. if the jacobian is not banded, set
//         mu to at least n - 1.
//
//       wa1 and wa2 are work arrays of length n. if ml + mu + 1 is at
//         least n, then the jacobian is considered dense, and wa2 is
//         not referenced.
{
  double eps;
  double epsmch;
  double h;
  int i;
  int j;
  int k;
  int msum;
  double temp;
//
//  EPSMCH is the machine precision.
//
  epsmch = std::numeric_limits<double>::epsilon();

  eps = sqrt ( fmax ( epsfcn, epsmch ) );
  msum = ml + mu + 1;
//
//  Computation of dense approximate jacobian.
//
  if ( n <= msum )
  {
    for ( j = 0; j < n; j++ )
    {
      temp = x[j];
      h = eps * fabs ( temp );
      if ( h == 0.0 )
      {
        h = eps;
      }
      x[j] = temp + h;
      fcn ( n, x, wa1 );
      x[j] = temp;
      for ( i = 0; i < n; i++ )
      {
        fjac[i+j*ldfjac] = ( wa1[i] - fvec[i] ) / h;
      }
    }
  }
//
//  Computation of a banded approximate jacobian.
//
  else
  {
    for ( k = 0; k < msum; k++ )
    {
      for ( j = k; j < n; j = j + msum )
      {
        wa2[j] = x[j];
        h = eps * fabs ( wa2[j] );
        if ( h == 0.0 )
        {
          h = eps;
        }
        x[j] = wa2[j] + h;
      }
      fcn ( n, x, wa1 );
      for ( j = k; j < n; j = j + msum )
      {
        x[j] = wa2[j];
        h = eps * fabs ( wa2[j] );
        if ( h == 0.0 )
        {
          h = eps;
        }
        for ( i = 0; i < n; i++ )
        {
          if ( j - mu <= i && i <= j + ml )
          {
            fjac[i+j*ldfjac] = ( wa1[i] - fvec[i] ) / h;
          }
          else
          {
            fjac[i+j*ldfjac] = 0.0;
          }
        }
      }
    }
  }
  return;
}
//****************************************************************************80

void fdjac_bdf2 ( void dydt ( double t, double x[], double f[] ),
  int n, double t1, double x1[], double t2, double x2[], double t3, double x3[], 
  double fvec[], double fjac[], int ldfjac, int ml, int mu, double epsfcn, 
  double wa1[], double wa2[] )

//****************************************************************************80
//
//  Purpose:
//
//    fdjac_bdf2() estimates a Jacobian matrix using forward differences.
//
//  Discussion:
//
//    This function computes a forward-difference approximation
//    to the N by N jacobian matrix associated with a specified
//    problem of N functions in N variables. 
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    18 November 2023
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    This version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Input:
//
//    void dydt ( double t, double y[], double dy[] ):
//    the name of the user-supplied code which
//    evaluates the right hand side of the ODE.
//
//    int n: the number of functions and variables.
//
//    double t1, x1[n], t2, x2[n], t3, x3[n]: three sets of
//    data at a sequence of times.
//
//    double fvec[n]: the functions evaluated at x3.
//
//    int ldfjac: the leading dimension of the array fjac, not less than N.
//
//    int ml, mu: the number of subdiagonals and
//    superdiagonals within the band of the jacobian matrix.  If the
//    jacobian is not banded, set ML and MU to N-1.
//
//    double epsfcn: determines a suitable step length for the 
//    forward-difference approximation. this approximation assumes 
//    that the relative errors in the functions are of the order 
//    of epsfcn.
//
//  Output:
//
//    double fjac[n*n]: the approximate jacobian matrix at x3.
//
//  Workspace:
//
//    double wa1[n], wa2[n]:  if ml + mu + 1 is at least n, then 
//    the jacobian is considered dense, and wa2 is not referenced.
//
{
  double eps;
  double epsmch;
  double h;
  int i;
  int j;
  int k;
  int msum;
  double temp;
//
//  EPSMCH is the machine precision.
//
  epsmch = std::numeric_limits<double>::epsilon();

  eps = sqrt ( fmax ( epsfcn, epsmch ) );
  msum = ml + mu + 1;
//
//  Computation of dense approximate jacobian.
//
  if ( n <= msum )
  {
    for ( j = 0; j < n; j++ )
    {
      temp = x3[j];
      h = eps * fabs ( temp );
      if ( h == 0.0 )
      {
        h = eps;
      }
      x3[j] = temp + h;
      bdf2_residual ( dydt, n, t1, x1, t2, x2, t3, x3, wa1 );
      x3[j] = temp;
      for ( i = 0; i < n; i++ )
      {
        fjac[i+j*ldfjac] = ( wa1[i] - fvec[i] ) / h;
      }
    }
  }
//
//  Computation of a banded approximate jacobian.
//
  else
  {
    for ( k = 0; k < msum; k++ )
    {
      for ( j = k; j < n; j = j + msum )
      {
        wa2[j] = x3[j];
        h = eps * fabs ( wa2[j] );
        if ( h == 0.0 )
        {
          h = eps;
        }
        x3[j] = wa2[j] + h;
      }
      bdf2_residual ( dydt, n, t1, x1, t2, x2, t3, x3, wa1 );
      for ( j = k; j < n; j = j + msum )
      {
        x3[j] = wa2[j];
        h = eps * fabs ( wa2[j] );
        if ( h == 0.0 )
        {
          h = eps;
        }
        for ( i = 0; i < n; i++ )
        {
          if ( j - mu <= i && i <= j + ml )
          {
            fjac[i+j*ldfjac] = ( wa1[i] - fvec[i] ) / h;
          }
          else
          {
            fjac[i+j*ldfjac] = 0.0;
          }
        }
      }
    }
  }
  return;
}
//****************************************************************************80

void fdjac_be ( void dydt ( double t, double x[], double f[] ),
  int n, double to, double xo[], double t, double x[], double fvec[], 
  double fjac[], int ldfjac, int ml, int mu, double epsfcn, double wa1[], 
  double wa2[] )

//****************************************************************************80
//
//  Purpose:
//
//    fdjac_be() estimates a Jacobian matrix using forward differences.
//
//  Discussion:
//
//    This function computes a forward-difference approximation
//    to the N by N jacobian matrix associated with a specified
//    problem of N functions in N variables. 
//
//    The original code fdjac1() was modified to deal with problems
//    involving a backward Euler residual.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    09 November 2023
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    This version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Input:
//
//    void dydt ( double t, double y[], double dy[] ):
//    the name of the user-supplied code which
//    evaluates the right hand side of the ODE.
//
//    int n: the number of functions and variables.
//
//    double to, xo[n]: the old time and solution.
//   
//    double t, x[n]: the new time and current solution estimate. 
//
//    double fvec[n]: the functions evaluated at X.
//
//    int ldfjac: the leading dimension of the array fjac, not less than N.
//
//    int ml, mu: the number of subdiagonals and
//    superdiagonals within the band of the jacobian matrix.  If the
//    jacobian is not banded, set ML and MU to N-1.
//
//    double epsfcn: determines a suitable step length for the 
//    forward-difference approximation. this approximation assumes 
//    that the relative errors in the functions are of the order 
//    of epsfcn.
//
//  Output:
//
//    double fjac[n*n]: the approximate jacobian matrix at X.
//
//  Workspace:
//
//    double wa1[n], wa2[n]:  if ml + mu + 1 is at least n, then 
//    the jacobian is considered dense, and wa2 is not referenced.
//
{
  double eps;
  double epsmch;
  double h;
  int i;
  int j;
  int k;
  int msum;
  double temp;
//
//  EPSMCH is the machine precision.
//
  epsmch = std::numeric_limits<double>::epsilon();

  eps = sqrt ( fmax ( epsfcn, epsmch ) );
  msum = ml + mu + 1;
//
//  Computation of dense approximate jacobian.
//
  if ( n <= msum )
  {
    for ( j = 0; j < n; j++ )
    {
      temp = x[j];
      h = eps * fabs ( temp );
      if ( h == 0.0 )
      {
        h = eps;
      }
      x[j] = temp + h;
      backward_euler_residual ( dydt, n, to, xo, t, x, wa1 );
      x[j] = temp;
      for ( i = 0; i < n; i++ )
      {
        fjac[i+j*ldfjac] = ( wa1[i] - fvec[i] ) / h;
      }
    }
  }
//
//  Computation of a banded approximate jacobian.
//
  else
  {
    for ( k = 0; k < msum; k++ )
    {
      for ( j = k; j < n; j = j + msum )
      {
        wa2[j] = x[j];
        h = eps * fabs ( wa2[j] );
        if ( h == 0.0 )
        {
          h = eps;
        }
        x[j] = wa2[j] + h;
      }
      backward_euler_residual ( dydt, n, to, xo, t, x, wa1 );
      for ( j = k; j < n; j = j + msum )
      {
        x[j] = wa2[j];
        h = eps * fabs ( wa2[j] );
        if ( h == 0.0 )
        {
          h = eps;
        }
        for ( i = 0; i < n; i++ )
        {
          if ( j - mu <= i && i <= j + ml )
          {
            fjac[i+j*ldfjac] = ( wa1[i] - fvec[i] ) / h;
          }
          else
          {
            fjac[i+j*ldfjac] = 0.0;
          }
        }
      }
    }
  }
  return;
}
//****************************************************************************80

void fdjac_tr ( void dydt ( double t, double x[], double f[] ),
  int n, double to, double xo[], double tn, double xn[], double fvec[], 
  double fjac[], int ldfjac, int ml, int mu, double epsfcn, double wa1[], 
  double wa2[] )

//****************************************************************************80
//
//  Purpose:
//
//    fdjac_tr() estimates a Jacobian matrix using forward differences.
//
//  Discussion:
//
//    This function computes a forward-difference approximation
//    to the N by N jacobian matrix associated with a specified
//    problem of N functions in N variables. 
//
//    The original code fdjac1() was modified to deal with problems
//    involving a trapezoidal ODE solver residual.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    16 November 2023
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    This version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Input:
//
//    void dydt ( double t, double y[], double dy[] ):
//    the name of the user-supplied code which
//    evaluates the right hand side of the ODE.
//
//    int n: the number of functions and variables.
//
//    double to, xo[n]: the old time and solution.
//   
//    double tn, xn[n]: the new time and current solution estimate. 
//
//    double fvec[n]: the functions evaluated at X.
//
//    int ldfjac: the leading dimension of the array fjac, not less than N.
//
//    int ml, mu: the number of subdiagonals and
//    superdiagonals within the band of the jacobian matrix.  If the
//    jacobian is not banded, set ML and MU to N-1.
//
//    double epsfcn: determines a suitable step length for the 
//    forward-difference approximation. this approximation assumes 
//    that the relative errors in the functions are of the order 
//    of epsfcn.
//
//  Output:
//
//    double fjac[n*n]: the approximate jacobian matrix at X.
//
//  Workspace:
//
//    double wa1[n], wa2[n]:  if ml + mu + 1 is at least n, then 
//    the jacobian is considered dense, and wa2 is not referenced.
//
{
  double eps;
  double epsmch;
  double h;
  int i;
  int j;
  int k;
  int msum;
  double temp;
//
//  EPSMCH is the machine precision.
//
  epsmch = std::numeric_limits<double>::epsilon();

  eps = sqrt ( fmax ( epsfcn, epsmch ) );
  msum = ml + mu + 1;
//
//  Computation of dense approximate jacobian.
//
  if ( n <= msum )
  {
    for ( j = 0; j < n; j++ )
    {
      temp = xn[j];
      h = eps * fabs ( temp );
      if ( h == 0.0 )
      {
        h = eps;
      }
      xn[j] = temp + h;
      trapezoidal_residual ( dydt, n, to, xo, tn, xn, wa1 );
      xn[j] = temp;
      for ( i = 0; i < n; i++ )
      {
        fjac[i+j*ldfjac] = ( wa1[i] - fvec[i] ) / h;
      }
    }
  }
//
//  Computation of a banded approximate jacobian.
//
  else
  {
    for ( k = 0; k < msum; k++ )
    {
      for ( j = k; j < n; j = j + msum )
      {
        wa2[j] = xn[j];
        h = eps * fabs ( wa2[j] );
        if ( h == 0.0 )
        {
          h = eps;
        }
        xn[j] = wa2[j] + h;
      }
      trapezoidal_residual ( dydt, n, to, xo, tn, xn, wa1 );
      for ( j = k; j < n; j = j + msum )
      {
        xn[j] = wa2[j];
        h = eps * fabs ( wa2[j] );
        if ( h == 0.0 )
        {
          h = eps;
        }
        for ( i = 0; i < n; i++ )
        {
          if ( j - mu <= i && i <= j + ml )
          {
            fjac[i+j*ldfjac] = ( wa1[i] - fvec[i] ) / h;
          }
          else
          {
            fjac[i+j*ldfjac] = 0.0;
          }
        }
      }
    }
  }
  return;
}
//****************************************************************************80

int fsolve ( void fcn ( int n, double x[], double fvec[] ), int n,
  double x[], double fvec[], double tol, double wa[], int lwa )

//****************************************************************************80
//
//  Purpose:
//
//    fsolve() finds a zero of a system of N nonlinear equations. 
//
//  Discussion:
//
//    The code finds a zero of a system of
//    N nonlinear functions in N variables by a modification
//    of the Powell hybrid method.  
//
//    This is done by using the more general nonlinear equation solver HYBRD.  
//
//    The user must provide FCN, which calculates the functions.
//
//    The jacobian is calculated by a forward-difference approximation.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    26 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    This version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Parameters:
//
//    Input, void FCN(int n,double x[], double fx[]): the name of the 
//    C++ routine which returns in fx[] the function value at the n-dimensional
//    vector x[].
//
//    Input, int N, the number of functions and variables.
//
//    Input/output, double X[N].  On input, an initial estimate of the solution.
//    On output, the final estimate of the solution.
//
//    Output, double FVEC[N], the functions evaluated at the output X.
//
//    Input, double TOL, a nonnegative variable. tTermination occurs when the 
//    algorithm estimates that the relative error between X and the solution 
//    is at most TOL.
//
//       info is an integer output variable.
//         info is set as follows.
//         info = 0   improper input parameters.
//         info = 1   algorithm estimates that the relative error
//                    between x and the solution is at most tol.
//         info = 2   number of calls to fcn has reached or exceeded
//                    200*(n+1).
//         info = 3   tol is too small. no further improvement in
//                    the approximate solution x is possible.
//         info = 4   iteration is not making good progress.
//
//    Workspace, double WA[LWA].
//
//       lwa is a positive integer input variable not less than
//         (n*(3*n+13))/2.
//
{
  double epsfcn;
  double factor;
  int index;
  int info;
  int j;
  int lr;
  int maxfev;
  int ml;
  int mode;
  int mu;
  int nfev;
  double xtol;

  info = 0;
//
//  Check the input.
//
  if ( n <= 0 )
  {
    return info;
  }
  if ( tol <= 0.0 )
  {
    return info;
  }
  if ( lwa < ( n * ( 3 * n + 13 ) ) / 2 )
  {
    return info;
  }
//
//  Call HYBRD.
//
  xtol = tol;
  maxfev = 200 * ( n + 1 );
  ml = n - 1;
  mu = n - 1;
  epsfcn = 0.0;
  for ( j = 0; j < n; j++ )
  {
    wa[j] = 1.0;
  }
  mode = 2;
  factor = 100.0;
  nfev = 0;
  lr = ( n * ( n + 1 ) ) / 2;
  index = 6 * n + lr;

  info = hybrd ( fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, wa, mode, 
    factor, nfev, wa+index, n, wa+6*n, lr, wa+n, wa+2*n, wa+3*n, wa+4*n, 
    wa+5*n );

  if ( info == 5 )
  {
    info = 4;
  }
  return info;
}
//****************************************************************************80

int fsolve_bdf2 ( void dydt ( double t, double x[], double f[] ), int n, 
  double t1, double x1[], double t2, double x2[], double t3, double x3[], 
  double fvec[], double tol )

//****************************************************************************80
//
//  Purpose:
//
//    fsolve_bdf2() finds a zero of a system of N nonlinear equations. 
//
//  Discussion:
//
//    The code finds a zero of a system of
//    N nonlinear functions in N variables by a modification
//    of the Powell hybrid method.  
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    18 November 2023
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    This version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Input:
//
//    void dydt ( double t, double y[], double dy[] ):
//    the name of the user-supplied code which
//    evaluates the right hand side of the ODE.
//
//    int n: the number of functions and variables.
//
//    double t1, x1[n], t2, x2[n], t3, x3[n]: three sets of
//    data at a sequence of times. 
//
//    double tol:  Satisfactory termination occurs when the algorithm
//    estimates that the relative error between X and the solution is at
//    most TOL.  TOL should be nonnegative.
//
//  Output:
//
//    double x3[n]: the final estimate of the solution.
//
//    double fvec[n]: the functions evaluated at the output X.
//
//    int info: status flag.  A value of 1 represents success.
//    0: improper input parameters.
//    1: algorithm estimates that the relative error
//       between x and the solution is at most tol.
//    2: number of calls to fcn has reached or exceeded 200*(n+1).
//    3: tol is too small. no further improvement in
//       the approximate solution x is possible.
//    4: iteration is not making good progress.
//
{
  double epsfcn;
  double factor;
  int index;
  int info;
  int j;
  int lr;
  int lwa;
  int maxfev;
  int ml;
  int mode;
  int mu;
  int nfev;
  double *wa;
  double xtol;

  info = 0;
//
//  Check the input.
//
  if ( n <= 0 )
  {
    return info;
  }
  if ( tol <= 0.0 )
  {
    return info;
  }

  lwa = ( n * ( 3 * n + 13 ) ) / 2;
  wa = new double[lwa];
//
//  Call HYBRD.
//
  xtol = tol;
  maxfev = 200 * ( n + 1 );
  ml = n - 1;
  mu = n - 1;
  epsfcn = 0.0;
  for ( j = 0; j < n; j++ )
  {
    wa[j] = 1.0;
  }
  mode = 2;
  factor = 100.0;
  nfev = 0;
  lr = ( n * ( n + 1 ) ) / 2;
  index = 6 * n + lr;

  info = hybrd_bdf2 ( dydt, n, t1, x1, t2, x2, t3, x3, fvec, xtol, maxfev, 
    ml, mu, epsfcn, wa, mode, factor, nfev, wa+index, n, wa+6*n, lr,
    wa+n, wa+2*n, wa+3*n, wa+4*n, wa+5*n );

  if ( info == 5 )
  {
    info = 4;
  }

  delete [] wa;

  return info;
}
//****************************************************************************80

int fsolve_be ( void dydt ( double t, double x[], double f[] ), int n, 
  double to, double xo[], double t, double x[], double fvec[], double tol )

//****************************************************************************80
//
//  Purpose:
//
//    fsolve_be() finds a zero of a system of N nonlinear equations. 
//
//  Discussion:
//
//    The code finds a zero of a system of
//    N nonlinear functions in N variables by a modification
//    of the Powell hybrid method.  
//
//    The original code fsolve() was modified to deal with problems
//    involving a backward Euler residual.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    09 November 2023
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    This version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Input:
//
//    void dydt ( double t, double y[], double dy[] ):
//    the name of the user-supplied code which
//    evaluates the right hand side of the ODE.
//
//    int n: the number of functions and variables.
//
//    double to, xo[n]: the old time and solution.
//   
//    double t, x[n]: the new time and current solution estimate. 
//
//    double tol:  Satisfactory termination occurs when the algorithm
//    estimates that the relative error between X and the solution is at
//    most TOL.  TOL should be nonnegative.
//
//  Output:
//
//    double x[n]: the final estimate of the solution.
//
//    double fvec[n]: the functions evaluated at the output X.
//
//    int info: status flag.  A value of 1 represents success.
//    0: improper input parameters.
//    1: algorithm estimates that the relative error
//       between x and the solution is at most tol.
//    2: number of calls to fcn has reached or exceeded 200*(n+1).
//    3: tol is too small. no further improvement in
//       the approximate solution x is possible.
//    4: iteration is not making good progress.
//
{
  double epsfcn;
  double factor;
  int index;
  int info;
  int j;
  int lr;
  int lwa;
  int maxfev;
  int ml;
  int mode;
  int mu;
  int nfev;
  double *wa;
  double xtol;

  info = 0;
//
//  Check the input.
//
  if ( n <= 0 )
  {
    return info;
  }
  if ( tol <= 0.0 )
  {
    return info;
  }

  lwa = ( n * ( 3 * n + 13 ) ) / 2;
  wa = new double[lwa];
//
//  Call HYBRD.
//
  xtol = tol;
  maxfev = 200 * ( n + 1 );
  ml = n - 1;
  mu = n - 1;
  epsfcn = 0.0;
  for ( j = 0; j < n; j++ )
  {
    wa[j] = 1.0;
  }
  mode = 2;
  factor = 100.0;
  nfev = 0;
  lr = ( n * ( n + 1 ) ) / 2;
  index = 6 * n + lr;

  info = hybrd_be ( dydt, n, to, xo, t, x, fvec, xtol, maxfev, ml, mu, 
    epsfcn, wa, mode, factor, nfev, wa+index, n, wa+6*n, lr,
    wa+n, wa+2*n, wa+3*n, wa+4*n, wa+5*n );

  if ( info == 5 )
  {
    info = 4;
  }

  delete [] wa;

  return info;
}
//****************************************************************************80

int fsolve_tr ( void dydt ( double t, double x[], double f[] ), int n, 
  double to, double xo[], double tn, double xn[], double fvec[], double tol )

//****************************************************************************80
//
//  Purpose:
//
//    fsolve_tr() finds a zero of a system of N nonlinear equations. 
//
//  Discussion:
//
//    The code finds a zero of a system of
//    N nonlinear functions in N variables by a modification
//    of the Powell hybrid method.  
//
//    The original code fsolve() was modified to deal with problems
//    involving a trapezoidal ODE solver residual.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    16 November 2023
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    This version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Input:
//
//    void dydt ( double t, double y[], double dy[] ):
//    the name of the user-supplied code which
//    evaluates the right hand side of the ODE.
//
//    int n: the number of functions and variables.
//
//    double to, xo[n]: the old time and solution.
//   
//    double tn, xn[n]: the new time and current solution estimate. 
//
//    double tol:  Satisfactory termination occurs when the algorithm
//    estimates that the relative error between X and the solution is at
//    most TOL.  TOL should be nonnegative.
//
//  Output:
//
//    double xn[n]: the final estimate of the solution.
//
//    double fvec[n]: the functions evaluated at the output xn.
//
//    int info: status flag.  A value of 1 represents success.
//    0: improper input parameters.
//    1: algorithm estimates that the relative error
//       between x and the solution is at most tol.
//    2: number of calls to fcn has reached or exceeded 200*(n+1).
//    3: tol is too small. no further improvement in
//       the approximate solution x is possible.
//    4: iteration is not making good progress.
//
{
  double epsfcn;
  double factor;
  int index;
  int info;
  int j;
  int lr;
  int lwa;
  int maxfev;
  int ml;
  int mode;
  int mu;
  int nfev;
  double *wa;
  double xtol;

  info = 0;
//
//  Check the input.
//
  if ( n <= 0 )
  {
    return info;
  }
  if ( tol <= 0.0 )
  {
    return info;
  }

  lwa = ( n * ( 3 * n + 13 ) ) / 2;
  wa = new double[lwa];
//
//  Call HYBRD.
//
  xtol = tol;
  maxfev = 200 * ( n + 1 );
  ml = n - 1;
  mu = n - 1;
  epsfcn = 0.0;
  for ( j = 0; j < n; j++ )
  {
    wa[j] = 1.0;
  }
  mode = 2;
  factor = 100.0;
  nfev = 0;
  lr = ( n * ( n + 1 ) ) / 2;
  index = 6 * n + lr;

  info = hybrd_tr ( dydt, n, to, xo, tn, xn, fvec, xtol, maxfev, ml, mu, 
    epsfcn, wa, mode, factor, nfev, wa+index, n, wa+6*n, lr,
    wa+n, wa+2*n, wa+3*n, wa+4*n, wa+5*n );

  if ( info == 5 )
  {
    info = 4;
  }

  delete [] wa;

  return info;
}
//****************************************************************************80

int hybrd ( void fcn ( int n, double x[], double fvec[] ),
  int n, double x[],
  double fvec[], double xtol, int maxfev, int ml, int mu, double epsfcn,
  double diag[], int mode, double factor, int nfev,
  double fjac[], int ldfjac, double r[], int lr, double qtf[], double wa1[],
  double wa2[], double wa3[], double wa4[] )

//****************************************************************************80
//
//  Purpose:
//
//    hybrd() finds a zero of a system of N nonlinear equations.
//
//  Discussion:
//
//    The code finds a zero of a system of
//    N nonlinear functions in N variables by a modification
//    of the Powell hybrid method. 
//
//    The user must provide FCN, which calculates the functions. 
//
//    The jacobian is calculated by a forward-difference approximation.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    08 April 2010
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    This version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Parameters:
//
//    Input, void FCN(int n,double x[], double fx[]): the name of the 
//    C++ routine which returns in fx[] the function value at the n-dimensional
//    vector x[].
//
//    Input, int N, the number of functions and variables.
//
//    Input/output, double X[N].  On input an initial estimate of the solution.
//    On output, the final estimate of the solution.
//
//    Output, double FVEC[N], the functions evaluated at the output value of X.
//
//    Input, double XTOL, a nonnegative value.  Termination occurs when the 
//    relative error between two consecutive iterates is at most XTOL.
//
//    Input, int MAXFEV.  Termination occurs when the number of calls to FCN 
//    is at least MAXFEV by the end of an iteration.
//
//    Input, int ML, specifies the number of subdiagonals within the band of 
//    the jacobian matrix.  If the jacobian is not banded, set
//    ml to at least n - 1.
//
//    Input, int MU, specifies the number of superdiagonals within the band 
//    of the jacobian matrix.  If the jacobian is not banded, set
//    mu to at least n - 1.
//
//       epsfcn is an input variable used in determining a suitable
//         step length for the forward-difference approximation. this
//         approximation assumes that the relative errors in the
//         functions are of the order of epsfcn. if epsfcn is less
//         than the machine precision, it is assumed that the relative
//         errors in the functions are of the order of the machine
//         precision.
//
//       diag is an array of length n. if mode = 1 (see
//         below), diag is internally set. if mode = 2, diag
//         must contain positive entries that serve as
//         multiplicative scale factors for the variables.
//
//       mode is an integer input variable. if mode = 1, the
//         variables will be scaled internally. if mode = 2,
//         the scaling is specified by the input diag. other
//         values of mode are equivalent to mode = 1.
//
//       factor is a positive input variable used in determining the
//         initial step bound. this bound is set to the product of
//         factor and the euclidean norm of diag*x if nonzero, or else
//         to factor itself. in most cases factor should lie in the
//         interval (.1,100.). 100. is a generally recommended value.
//
//       info is an integer output variable.
//         info = 0   improper input parameters.
//         info = 1   relative error between two consecutive iterates
//                    is at most xtol.
//         info = 2   number of calls to fcn has reached or exceeded
//                    maxfev.
//         info = 3   xtol is too small. no further improvement in
//                    the approximate solution x is possible.
//         info = 4   iteration is not making good progress, as
//                    measured by the improvement from the last
//                    five jacobian evaluations.
//         info = 5   iteration is not making good progress, as
//                    measured by the improvement from the last
//                    ten iterations.
//
//       nfev is an integer output variable set to the number of
//         calls to fcn.
//
//       fjac is an output n by n array which contains the
//         orthogonal matrix q produced by the qr factorization
//         of the final approximate jacobian.
//
//       ldfjac is a positive integer input variable not less than n
//         which specifies the leading dimension of the array fjac.
//
//       r is an output array of length lr which contains the
//         upper triangular matrix produced by the qr factorization
//         of the final approximate jacobian, stored rowwise.
//
//       lr is a positive integer input variable not less than
//         (n*(n+1))/2.
//
//       qtf is an output array of length n which contains
//         the vector (q transpose)*fvec.
//
//       wa1, wa2, wa3, and wa4 are work arrays of length n.
//
{
  double actred;
  double delta;
  double epsmch;
  double fnorm;
  double fnorm1;
  int i;
  int info;
  int iter;
  int iwa[1];
  int j;
  bool jeval;
  int l;
  int msum;
  int ncfail;
  int ncsuc;
  int nslow1;
  int nslow2;
  const double p001 = 0.001;
  const double p0001 = 0.0001;
  const double p1 = 0.1;
  const double p5 = 0.5;
  double pnorm;
  double prered;
  double ratio;
  double sum;
  double temp;
  double xnorm;
//
//  Certain loops in this function were kept closer to their original FORTRAN77
//  format, to avoid confusing issues with the array index L.  These loops are
//  marked "DO NOT ADJUST", although they certainly could be adjusted (carefully)
//  once the initial translated code is tested.
//

//
//  EPSMCH is the machine precision.
//
  epsmch = std::numeric_limits<double>::epsilon();

  info = 0;
  nfev = 0;
//
//  Check the input parameters.
//
  if ( n <= 0 )
  {
    info = 0;
    return info;
  }
  if ( xtol < 0.0 )
  {
    info = 0;
    return info;
  }
  if ( maxfev <= 0 )
  {
    info = 0;
    return info;
  }
  if ( ml < 0 )
  {
    info = 0;
    return info;
  }
  if ( mu < 0 )
  {
    info = 0;
    return info;
  }
  if ( factor <= 0.0 )
  {
    info = 0;
    return info;
  }
  if ( ldfjac < n )
  {
    info = 0;
    return info;
  }
  if ( lr < ( n * ( n + 1 ) ) / 2 )
  {
    info = 0;
    return info;
  }
  if ( mode == 2 )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( diag[j] <= 0.0 )
      {
        info = 0;
        return info;
      }
    }
  }
//
//  Evaluate the function at the starting point and calculate its norm.
//
  fcn ( n, x, fvec );
  nfev = 1;

  fnorm = enorm ( n, fvec );
//
//  Determine the number of calls to FCN needed to compute the jacobian matrix.
//
  if ( ml + mu + 1 < n )
  {
    msum = ml + mu + 1;
  }
  else
  {
    msum = n;
  }
//
//  Initialize iteration counter and monitors.
//
  iter = 1;
  ncsuc = 0;
  ncfail = 0;
  nslow1 = 0;
  nslow2 = 0;
//
//  Beginning of the outer loop.
//
  for ( ; ; )
  {
    jeval = true;
//
//  Calculate the jacobian matrix.
//
    fdjac1 ( fcn, n, x, fvec, fjac, ldfjac, ml, mu, epsfcn, wa1, wa2 );

    nfev = nfev + msum;
//
//  Compute the QR factorization of the jacobian.
//
    qrfac ( n, n, fjac, ldfjac, false, iwa, 1, wa1, wa2 );
//
//  On the first iteration and if MODE is 1, scale according
//  to the norms of the columns of the initial jacobian.
//
    if ( iter == 1 )
    {
      if ( mode == 1 )
      {
        for ( j = 0; j < n; j++ )
        {
          if ( wa2[j] != 0.0 )
          {
            diag[j] = wa2[j];
          }
          else
          {
            diag[j] = 1.0;
          }
        }
      }
//
//  On the first iteration, calculate the norm of the scaled X
//  and initialize the step bound DELTA.
//
      for ( j = 0; j < n; j++ )
      {
        wa3[j] = diag[j] * x[j];
      }
      xnorm = enorm ( n, wa3 );

      if ( xnorm == 0.0 )
      {
        delta = factor;
      }
      else
      {
        delta = factor * xnorm;
      }
    }
//
//  Form Q' * FVEC and store in QTF.
//
    for ( i = 0; i < n; i++ )
    {
      qtf[i] = fvec[i];
    }
    for ( j = 0; j < n; j++ )
    {
      if ( fjac[j+j*ldfjac] != 0.0 )
      {
        sum = 0.0;
        for ( i = j; i < n; i++ )
        {
          sum = sum + fjac[i+j*ldfjac] * qtf[i];
        }
        temp = - sum / fjac[j+j*ldfjac];
        for ( i = j; i < n; i++ )
        {
          qtf[i] = qtf[i] + fjac[i+j*ldfjac] * temp;
        }
      }
    }
//
//  Copy the triangular factor of the QR factorization into R.
//
//  DO NOT ADJUST THIS LOOP, BECAUSE OF L.
//
    for ( j = 1; j <= n; j++ )
    {
      l = j;
      for ( i = 1; i <= j - 1; i++ )
      {
        r[l-1] = fjac[(i-1)+(j-1)*ldfjac];
        l = l + n - i;
      }
      r[l-1] = wa1[j-1];
      if ( wa1[j-1] == 0.0 )
      {
        cerr << "  Matrix is singular.\n";
      }
    }
//
//  Accumulate the orthogonal factor in FJAC.
//
    qform ( n, n, fjac, ldfjac );
//
//  Rescale if necessary.
//
    if ( mode == 1 )
    {
      for ( j = 0; j < n; j++ )
      {
        diag[j] = fmax ( diag[j], wa2[j] );
      }
    }
//
//  Beginning of the inner loop.
//
    for ( ; ; )
    {
//
//  Determine the direction P.
//
      dogleg ( n, r, lr, diag, qtf, delta, wa1, wa2, wa3 );
//
//  Store the direction P and X + P.  Calculate the norm of P.
//
      for ( j = 0; j < n; j++ )
      {
        wa1[j] = - wa1[j];
        wa2[j] = x[j] + wa1[j];
        wa3[j] = diag[j] * wa1[j];
      }
      pnorm = enorm ( n, wa3 );
//
//  On the first iteration, adjust the initial step bound.
//
      if ( iter == 1 )
      {
        delta = fmin ( delta, pnorm );
      }
//
//  Evaluate the function at X + P and calculate its norm.
//
      fcn ( n, wa2, wa4 );
      nfev = nfev + 1;
      fnorm1 = enorm ( n, wa4 );
//
//  Compute the scaled actual reduction.
//
      if ( fnorm1 < fnorm )
      {
        actred = 1.0 - ( fnorm1 / fnorm ) * ( fnorm1 / fnorm );
      }
      else
      {
        actred = - 1.0;
      }
//
//  Compute the scaled predicted reduction.
//
//  DO NOT ADJUST THIS LOOP, BECAUSE OF L.
//
      l = 1;
      for ( i = 1; i <= n; i++ )
      {
        sum = 0.0;
        for ( j = i; j <= n; j++ )
        {
          sum = sum + r[l-1] * wa1[j-1];
          l = l + 1;
        }
        wa3[i-1] = qtf[i-1] + sum;
      }
      temp = enorm ( n, wa3 );

      if ( temp < fnorm )
      {
        prered = 1.0 - ( temp / fnorm ) * ( temp / fnorm );
      }
      else
      {
        prered = 0.0;
      }
//
//  Compute the ratio of the actual to the predicted reduction.
//
      if ( 0.0 < prered )
      {
        ratio = actred / prered;
      }
      else
      {
        ratio = 0.0;
      }
//
//  Update the step bound.
//
      if ( ratio < p1 )
      {
        ncsuc = 0;
        ncfail = ncfail + 1;
        delta = p5 * delta;
      }
      else
      {
        ncfail = 0;
        ncsuc = ncsuc + 1;
        if ( p5 <= ratio || 1 < ncsuc )
        {
          delta = fmax ( delta, pnorm / p5 );
        }
        if ( fabs ( ratio - 1.0 ) <= p1 )
        {
          delta = pnorm / p5;
        }
      }
//
//  On successful iteration, update X, FVEC, and their norms.
//
      if ( p0001 <= ratio )
      {
        for ( j = 0; j < n; j++ )
        {
          x[j] = wa2[j];
          wa2[j] = diag[j] * x[j];
          fvec[j] = wa4[j];
        }
        xnorm = enorm ( n, wa2 );
        fnorm = fnorm1;
        iter = iter + 1;
      }
//
//  Determine the progress of the iteration.
//
      nslow1 = nslow1 + 1;
      if ( p001 <= actred )
      {
        nslow1 = 0;
      }
      if ( jeval )
      {
        nslow2 = nslow2 + 1;
      }
      if ( p1 <= actred )
      {
        nslow2 = 0;
      }
//
//  Test for convergence.
//
      if ( delta <= xtol * xnorm || fnorm == 0.0 )
      {
        info = 1;
        return info;
      }
//
//  Tests for termination and stringent tolerances.
//
      if ( maxfev <= nfev )
      {
        info = 2;
        return info;
      }
      if ( p1 * fmax ( p1 * delta, pnorm ) <= epsmch * xnorm )
      {
        info = 3;
        return info;
      }
      if ( nslow2 == 5 )
      {
        info = 4;
        return info;
      }
      if ( nslow1 == 10 )
      {
        info = 5;
        return info;
      }
//
//  Criterion for recalculating jacobian approximation by forward differences.
//
      if ( ncfail == 2 )
      {
        break;
      }
//
//  Calculate the rank one modification to the jacobian
//  and update QTF if necessary.
//
      for ( j = 0; j < n; j++ )
      {
        sum = 0.0;
        for ( i = 0; i < n; i++ )
        {
          sum = sum + fjac[i+j*ldfjac] * wa4[i];
        }
        wa2[j] = ( sum - wa3[j] ) / pnorm;
        wa1[j] = diag[j] * ( ( diag[j] * wa1[j] ) / pnorm );
        if ( p0001 <= ratio )
        {
          qtf[j] = sum;
        }
      }
//
//  Compute the QR factorization of the updated jacobian.
//
      r1updt ( n, n, r, lr, wa1, wa2, wa3 );
      r1mpyq ( n, n, fjac, ldfjac, wa2, wa3 );
      r1mpyq ( 1, n, qtf, 1, wa2, wa3 );

      jeval = false;
    }
//
//  End of the inner loop.
//
  }
//
//  End of the outer loop.
//
}
//****************************************************************************80

int hybrd_bdf2 ( void dydt ( double t, double x[], double f[] ), int n, 
  double t1, double x1[], double t2, double x2[], double t3, double x3[],
  double fvec[], double xtol, int maxfev, int ml, int mu, double epsfcn,
  double diag[], int mode, double factor, int nfev,
  double fjac[], int ldfjac, double r[], int lr, double qtf[], double wa1[],
  double wa2[], double wa3[], double wa4[] )

//****************************************************************************80
//
//  Purpose:
//
//    hybrd_bdf2() finds a zero of a system of N nonlinear equations.
//
//  Discussion:
//
//    The code finds a zero of a system of
//    N nonlinear functions in N variables by a modification
//    of the Powell hybrid method. 
//
//    The jacobian is calculated by a forward-difference approximation.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    18 November 2023
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    This version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Input:
//
//    void dydt ( double t, double y[], double dy[] ):
//    the name of the user-supplied code which
//    evaluates the right hand side of the ODE.
//
//    int n: the number of functions and variables.
//
//    double t1, x1[n], t2, x2[n], t3, x3[n]: three sets of
//    data at a sequence of times. 
//
//    double xtol: Termination occurs when the 
//    relative error between two consecutive iterates is at most XTOL.
//
//    int maxfev:  Termination occurs when the number of calls to FCN 
//    is at least MAXFEV by the end of an iteration.
//
//    int ml, mu: the number of subdiagonals and superdiagonals
//    within the band of the jacobian matrix.  If the jacobian is not 
//    banded, set ml and mu to at least n - 1.
//
//    double epsfcn: determines a suitable step length for the 
//    forward-difference approximation.  This approximation assumes that 
//    the relative errors in the functions are of the order of epsfcn. 
//    If epsfcn is less than the machine precision, it is assumed that 
//    the relative errors in the functions are of the order of the machine
//    precision.
//
//    double diag[n]: multiplicative scale factors for the 
//    variables.  Only needed as input if MODE = 2.
//
//    int mode: scaling option.
//    if mode = 1, the variables will be scaled internally. 
//    if mode = 2, the scaling is specified by the input diag.
//
//    double factor: determines the initial step bound. this bound is set 
//    to the product of factor and the euclidean norm of diag*x if nonzero, 
//    or else to factor itself. in most cases factor should lie in the
//    interval (.1,100.). 100. is a generally recommended value.
//
//    int ldfjac: not less than n, specifies the leading 
//    dimension of the array fjac.
//
//    int lr: the size of the R array, not less than (n*(n+1))/2.
//
//  Output:
//
//    double x[n]: the final estimate of the solution.
//
//    double fvec[n]: the functions evaluated at the output value of X.
//
//    double diag[n]: multiplicative scale factors for the 
//    variables.  Set internally if MODE = 1.
//
//    int INFO: status flag.  A value of 1 indicates success. 
//    0: improper input parameters.
//    1: algorithm estimates that the relative error
//       between x and the solution is at most tol.
//    2: number of calls to the deriviative has reached or exceeded 200*(n+1).
//    3: tol is too small. no further improvement in
//       the approximate solution x is possible.
//    4: iteration is not making good progress.
//
//    int nfev: the number of calls to the derivative function.
//
//    double fjac[n*n]: the orthogonal matrix q produced by the qr 
//    factorization of the final approximate jacobian.
//
//    double r[lr]: the upper triangular matrix produced by the qr 
//    factorization of the final approximate jacobian, stored rowwise.
//
//    double qtf[n]: contains the vector (q transpose)*fvec.
//
//  Workspace:
//
//    double wa1[n], wa2[n], wa3[n], wa4[n].
//
{
  double actred;
  double delta;
  double epsmch;
  double fnorm;
  double fnorm1;
  int i;
  int info;
  int iter;
  int iwa[1];
  int j;
  bool jeval;
  int l;
  int msum;
  int ncfail;
  int ncsuc;
  int nslow1;
  int nslow2;
  const double p001 = 0.001;
  const double p0001 = 0.0001;
  const double p1 = 0.1;
  const double p5 = 0.5;
  double pnorm;
  double prered;
  double ratio;
  double sum;
  double temp;
  double xnorm;
//
//  Certain loops in this function were kept closer to their original FORTRAN77
//  format, to avoid confusing issues with the array index L.  These loops are
//  marked "DO NOT ADJUST", although they certainly could be adjusted (carefully)
//  once the initial translated code is tested.
//

//
//  EPSMCH is the machine precision.
//
  epsmch = std::numeric_limits<double>::epsilon();

  info = 0;
  nfev = 0;
//
//  Check the input parameters.
//
  if ( n <= 0 )
  {
    info = 0;
    return info;
  }
  if ( xtol < 0.0 )
  {
    info = 0;
    return info;
  }
  if ( maxfev <= 0 )
  {
    info = 0;
    return info;
  }
  if ( ml < 0 )
  {
    info = 0;
    return info;
  }
  if ( mu < 0 )
  {
    info = 0;
    return info;
  }
  if ( factor <= 0.0 )
  {
    info = 0;
    return info;
  }
  if ( ldfjac < n )
  {
    info = 0;
    return info;
  }
  if ( lr < ( n * ( n + 1 ) ) / 2 )
  {
    info = 0;
    return info;
  }
  if ( mode == 2 )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( diag[j] <= 0.0 )
      {
        info = 0;
        return info;
      }
    }
  }
//
//  Evaluate the function at the starting point and calculate its norm.
//
  bdf2_residual ( dydt, n, t1, x1, t2, x2, t3, x3, fvec );
  nfev = 1;

  fnorm = enorm ( n, fvec );
//
//  Determine the number of calls to FCN needed to compute the jacobian matrix.
//
  if ( ml + mu + 1 < n )
  {
    msum = ml + mu + 1;
  }
  else
  {
    msum = n;
  }
//
//  Initialize iteration counter and monitors.
//
  iter = 1;
  ncsuc = 0;
  ncfail = 0;
  nslow1 = 0;
  nslow2 = 0;
//
//  Beginning of the outer loop.
//
  for ( ; ; )
  {
    jeval = true;
//
//  Calculate the jacobian matrix.
//
    fdjac_bdf2 ( dydt, n, t1, x1, t2, x2, t3, x3, fvec, fjac, ldfjac, 
      ml, mu, epsfcn, wa1, wa2 );

    nfev = nfev + msum;
//
//  Compute the QR factorization of the jacobian.
//
    qrfac ( n, n, fjac, ldfjac, false, iwa, 1, wa1, wa2 );
//
//  On the first iteration and if MODE is 1, scale according
//  to the norms of the columns of the initial jacobian.
//
    if ( iter == 1 )
    {
      if ( mode == 1 )
      {
        for ( j = 0; j < n; j++ )
        {
          if ( wa2[j] != 0.0 )
          {
            diag[j] = wa2[j];
          }
          else
          {
            diag[j] = 1.0;
          }
        }
      }
//
//  On the first iteration, calculate the norm of the scaled X
//  and initialize the step bound DELTA.
//
      for ( j = 0; j < n; j++ )
      {
        wa3[j] = diag[j] * x3[j];
      }
      xnorm = enorm ( n, wa3 );

      if ( xnorm == 0.0 )
      {
        delta = factor;
      }
      else
      {
        delta = factor * xnorm;
      }
    }
//
//  Form Q' * FVEC and store in QTF.
//
    for ( i = 0; i < n; i++ )
    {
      qtf[i] = fvec[i];
    }
    for ( j = 0; j < n; j++ )
    {
      if ( fjac[j+j*ldfjac] != 0.0 )
      {
        sum = 0.0;
        for ( i = j; i < n; i++ )
        {
          sum = sum + fjac[i+j*ldfjac] * qtf[i];
        }
        temp = - sum / fjac[j+j*ldfjac];
        for ( i = j; i < n; i++ )
        {
          qtf[i] = qtf[i] + fjac[i+j*ldfjac] * temp;
        }
      }
    }
//
//  Copy the triangular factor of the QR factorization into R.
//
//  DO NOT ADJUST THIS LOOP, BECAUSE OF L.
//
    for ( j = 1; j <= n; j++ )
    {
      l = j;
      for ( i = 1; i <= j - 1; i++ )
      {
        r[l-1] = fjac[(i-1)+(j-1)*ldfjac];
        l = l + n - i;
      }
      r[l-1] = wa1[j-1];
      if ( wa1[j-1] == 0.0 )
      {
        cerr << "hybrd_bdf2(): Matrix is singular.\n";
      }
    }
//
//  Accumulate the orthogonal factor in FJAC.
//
    qform ( n, n, fjac, ldfjac );
//
//  Rescale if necessary.
//
    if ( mode == 1 )
    {
      for ( j = 0; j < n; j++ )
      {
        diag[j] = fmax ( diag[j], wa2[j] );
      }
    }
//
//  Beginning of the inner loop.
//
    for ( ; ; )
    {
//
//  Determine the direction P.
//
      dogleg ( n, r, lr, diag, qtf, delta, wa1, wa2, wa3 );
//
//  Store the direction P and X + P.  Calculate the norm of P.
//
      for ( j = 0; j < n; j++ )
      {
        wa1[j] = - wa1[j];
        wa2[j] = x3[j] + wa1[j];
        wa3[j] = diag[j] * wa1[j];
      }
      pnorm = enorm ( n, wa3 );
//
//  On the first iteration, adjust the initial step bound.
//
      if ( iter == 1 )
      {
        delta = fmin ( delta, pnorm );
      }
//
//  Evaluate the function at X + P and calculate its norm.
//
      bdf2_residual ( dydt, n, t1, x1, t2, x2, t3, wa2, wa4 );
      nfev = nfev + 1;
      fnorm1 = enorm ( n, wa4 );
//
//  Compute the scaled actual reduction.
//
      if ( fnorm1 < fnorm )
      {
        actred = 1.0 - ( fnorm1 / fnorm ) * ( fnorm1 / fnorm );
      }
      else
      {
        actred = - 1.0;
      }
//
//  Compute the scaled predicted reduction.
//
//  DO NOT ADJUST THIS LOOP, BECAUSE OF L.
//
      l = 1;
      for ( i = 1; i <= n; i++ )
      {
        sum = 0.0;
        for ( j = i; j <= n; j++ )
        {
          sum = sum + r[l-1] * wa1[j-1];
          l = l + 1;
        }
        wa3[i-1] = qtf[i-1] + sum;
      }
      temp = enorm ( n, wa3 );

      if ( temp < fnorm )
      {
        prered = 1.0 - ( temp / fnorm ) * ( temp / fnorm );
      }
      else
      {
        prered = 0.0;
      }
//
//  Compute the ratio of the actual to the predicted reduction.
//
      if ( 0.0 < prered )
      {
        ratio = actred / prered;
      }
      else
      {
        ratio = 0.0;
      }
//
//  Update the step bound.
//
      if ( ratio < p1 )
      {
        ncsuc = 0;
        ncfail = ncfail + 1;
        delta = p5 * delta;
      }
      else
      {
        ncfail = 0;
        ncsuc = ncsuc + 1;
        if ( p5 <= ratio || 1 < ncsuc )
        {
          delta = fmax ( delta, pnorm / p5 );
        }
        if ( fabs ( ratio - 1.0 ) <= p1 )
        {
          delta = pnorm / p5;
        }
      }
//
//  On successful iteration, update X, FVEC, and their norms.
//
      if ( p0001 <= ratio )
      {
        for ( j = 0; j < n; j++ )
        {
          x3[j] = wa2[j];
          wa2[j] = diag[j] * x3[j];
          fvec[j] = wa4[j];
        }
        xnorm = enorm ( n, wa2 );
        fnorm = fnorm1;
        iter = iter + 1;
      }
//
//  Determine the progress of the iteration.
//
      nslow1 = nslow1 + 1;
      if ( p001 <= actred )
      {
        nslow1 = 0;
      }
      if ( jeval )
      {
        nslow2 = nslow2 + 1;
      }
      if ( p1 <= actred )
      {
        nslow2 = 0;
      }
//
//  Test for convergence.
//
      if ( delta <= xtol * xnorm || fnorm == 0.0 )
      {
        info = 1;
        return info;
      }
//
//  Tests for termination and stringent tolerances.
//
      if ( maxfev <= nfev )
      {
        info = 2;
        return info;
      }
      if ( p1 * fmax ( p1 * delta, pnorm ) <= epsmch * xnorm )
      {
        info = 3;
        return info;
      }
      if ( nslow2 == 5 )
      {
        info = 4;
        return info;
      }
      if ( nslow1 == 10 )
      {
        info = 5;
        return info;
      }
//
//  Criterion for recalculating jacobian approximation by forward differences.
//
      if ( ncfail == 2 )
      {
        break;
      }
//
//  Calculate the rank one modification to the jacobian
//  and update QTF if necessary.
//
      for ( j = 0; j < n; j++ )
      {
        sum = 0.0;
        for ( i = 0; i < n; i++ )
        {
          sum = sum + fjac[i+j*ldfjac] * wa4[i];
        }
        wa2[j] = ( sum - wa3[j] ) / pnorm;
        wa1[j] = diag[j] * ( ( diag[j] * wa1[j] ) / pnorm );
        if ( p0001 <= ratio )
        {
          qtf[j] = sum;
        }
      }
//
//  Compute the QR factorization of the updated jacobian.
//
      r1updt ( n, n, r, lr, wa1, wa2, wa3 );
      r1mpyq ( n, n, fjac, ldfjac, wa2, wa3 );
      r1mpyq ( 1, n, qtf, 1, wa2, wa3 );

      jeval = false;
    }
//
//  End of the inner loop.
//
  }
//
//  End of the outer loop.
//
}
//****************************************************************************80

int hybrd_be ( void dydt ( double t, double x[], double f[] ), int n, double to, 
  double xo[], double t, double x[],
  double fvec[], double xtol, int maxfev, int ml, int mu, double epsfcn,
  double diag[], int mode, double factor, int nfev,
  double fjac[], int ldfjac, double r[], int lr, double qtf[], double wa1[],
  double wa2[], double wa3[], double wa4[] )

//****************************************************************************80
//
//  Purpose:
//
//    hybrd_be() finds a zero of a system of N nonlinear equations.
//
//  Discussion:
//
//    The code finds a zero of a system of
//    N nonlinear functions in N variables by a modification
//    of the Powell hybrid method. 
//
//    The jacobian is calculated by a forward-difference approximation.
//
//    The original code hybrd() was modified to deal with problems
//    involving a backward Euler residual.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    09 November 2023
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    This version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Input:
//
//    void dydt ( double t, double y[], double dy[] ):
//    the name of the user-supplied code which
//    evaluates the right hand side of the ODE.
//
//    int n: the number of functions and variables.
//
//    double to, xo[n]: the old time and solution.
// 
//    double t, x[n]: the new time and current solution estimate. 
//
//    double xtol: Termination occurs when the 
//    relative error between two consecutive iterates is at most XTOL.
//
//    int maxfev:  Termination occurs when the number of calls to FCN 
//    is at least MAXFEV by the end of an iteration.
//
//    int ml, mu: the number of subdiagonals and superdiagonals
//    within the band of the jacobian matrix.  If the jacobian is not 
//    banded, set ml and mu to at least n - 1.
//
//    double epsfcn: determines a suitable step length for the 
//    forward-difference approximation.  This approximation assumes that 
//    the relative errors in the functions are of the order of epsfcn. 
//    If epsfcn is less than the machine precision, it is assumed that 
//    the relative errors in the functions are of the order of the machine
//    precision.
//
//    double diag[n]: multiplicative scale factors for the 
//    variables.  Only needed as input if MODE = 2.
//
//    int mode: scaling option.
//    if mode = 1, the variables will be scaled internally. 
//    if mode = 2, the scaling is specified by the input diag.
//
//    double factor: determines the initial step bound. this bound is set 
//    to the product of factor and the euclidean norm of diag*x if nonzero, 
//    or else to factor itself. in most cases factor should lie in the
//    interval (.1,100.). 100. is a generally recommended value.
//
//    int ldfjac: not less than n, specifies the leading 
//    dimension of the array fjac.
//
//    int lr: the size of the R array, not less than (n*(n+1))/2.
//
//  Output:
//
//    double x[n]: the final estimate of the solution.
//
//    double fvec[n]: the functions evaluated at the output value of X.
//
//    double diag[n]: multiplicative scale factors for the 
//    variables.  Set internally if MODE = 1.
//
//    int INFO: status flag.  A value of 1 indicates success. 
//    0: improper input parameters.
//    1: algorithm estimates that the relative error
//       between x and the solution is at most tol.
//    2: number of calls to the deriviative has reached or exceeded 200*(n+1).
//    3: tol is too small. no further improvement in
//       the approximate solution x is possible.
//    4: iteration is not making good progress.
//
//    int nfev: the number of calls to the derivative function.
//
//    double fjac[n*n]: the orthogonal matrix q produced by the qr 
//    factorization of the final approximate jacobian.
//
//    double r[lr]: the upper triangular matrix produced by the qr 
//    factorization of the final approximate jacobian, stored rowwise.
//
//    double qtf[n]: contains the vector (q transpose)*fvec.
//
//  Workspace:
//
//    double wa1[n], wa2[n], wa3[n], wa4[n].
//
{
  double actred;
  double delta;
  double epsmch;
  double fnorm;
  double fnorm1;
  int i;
  int info;
  int iter;
  int iwa[1];
  int j;
  bool jeval;
  int l;
  int msum;
  int ncfail;
  int ncsuc;
  int nslow1;
  int nslow2;
  const double p001 = 0.001;
  const double p0001 = 0.0001;
  const double p1 = 0.1;
  const double p5 = 0.5;
  double pnorm;
  double prered;
  double ratio;
  double sum;
  double temp;
  double xnorm;
//
//  Certain loops in this function were kept closer to their original FORTRAN77
//  format, to avoid confusing issues with the array index L.  These loops are
//  marked "DO NOT ADJUST", although they certainly could be adjusted (carefully)
//  once the initial translated code is tested.
//

//
//  EPSMCH is the machine precision.
//
  epsmch = std::numeric_limits<double>::epsilon();

  info = 0;
  nfev = 0;
//
//  Check the input parameters.
//
  if ( n <= 0 )
  {
    info = 0;
    return info;
  }
  if ( xtol < 0.0 )
  {
    info = 0;
    return info;
  }
  if ( maxfev <= 0 )
  {
    info = 0;
    return info;
  }
  if ( ml < 0 )
  {
    info = 0;
    return info;
  }
  if ( mu < 0 )
  {
    info = 0;
    return info;
  }
  if ( factor <= 0.0 )
  {
    info = 0;
    return info;
  }
  if ( ldfjac < n )
  {
    info = 0;
    return info;
  }
  if ( lr < ( n * ( n + 1 ) ) / 2 )
  {
    info = 0;
    return info;
  }
  if ( mode == 2 )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( diag[j] <= 0.0 )
      {
        info = 0;
        return info;
      }
    }
  }
//
//  Evaluate the function at the starting point and calculate its norm.
//
  backward_euler_residual ( dydt, n, to, xo, t, x, fvec );
  nfev = 1;

  fnorm = enorm ( n, fvec );
//
//  Determine the number of calls to FCN needed to compute the jacobian matrix.
//
  if ( ml + mu + 1 < n )
  {
    msum = ml + mu + 1;
  }
  else
  {
    msum = n;
  }
//
//  Initialize iteration counter and monitors.
//
  iter = 1;
  ncsuc = 0;
  ncfail = 0;
  nslow1 = 0;
  nslow2 = 0;
//
//  Beginning of the outer loop.
//
  for ( ; ; )
  {
    jeval = true;
//
//  Calculate the jacobian matrix.
//
    fdjac_be ( dydt, n, to, xo, t, x, fvec, fjac, ldfjac, ml, mu, epsfcn, 
      wa1, wa2 );

    nfev = nfev + msum;
//
//  Compute the QR factorization of the jacobian.
//
    qrfac ( n, n, fjac, ldfjac, false, iwa, 1, wa1, wa2 );
//
//  On the first iteration and if MODE is 1, scale according
//  to the norms of the columns of the initial jacobian.
//
    if ( iter == 1 )
    {
      if ( mode == 1 )
      {
        for ( j = 0; j < n; j++ )
        {
          if ( wa2[j] != 0.0 )
          {
            diag[j] = wa2[j];
          }
          else
          {
            diag[j] = 1.0;
          }
        }
      }
//
//  On the first iteration, calculate the norm of the scaled X
//  and initialize the step bound DELTA.
//
      for ( j = 0; j < n; j++ )
      {
        wa3[j] = diag[j] * x[j];
      }
      xnorm = enorm ( n, wa3 );

      if ( xnorm == 0.0 )
      {
        delta = factor;
      }
      else
      {
        delta = factor * xnorm;
      }
    }
//
//  Form Q' * FVEC and store in QTF.
//
    for ( i = 0; i < n; i++ )
    {
      qtf[i] = fvec[i];
    }
    for ( j = 0; j < n; j++ )
    {
      if ( fjac[j+j*ldfjac] != 0.0 )
      {
        sum = 0.0;
        for ( i = j; i < n; i++ )
        {
          sum = sum + fjac[i+j*ldfjac] * qtf[i];
        }
        temp = - sum / fjac[j+j*ldfjac];
        for ( i = j; i < n; i++ )
        {
          qtf[i] = qtf[i] + fjac[i+j*ldfjac] * temp;
        }
      }
    }
//
//  Copy the triangular factor of the QR factorization into R.
//
//  DO NOT ADJUST THIS LOOP, BECAUSE OF L.
//
    for ( j = 1; j <= n; j++ )
    {
      l = j;
      for ( i = 1; i <= j - 1; i++ )
      {
        r[l-1] = fjac[(i-1)+(j-1)*ldfjac];
        l = l + n - i;
      }
      r[l-1] = wa1[j-1];
      if ( wa1[j-1] == 0.0 )
      {
        cerr << "  Matrix is singular.\n";
      }
    }
//
//  Accumulate the orthogonal factor in FJAC.
//
    qform ( n, n, fjac, ldfjac );
//
//  Rescale if necessary.
//
    if ( mode == 1 )
    {
      for ( j = 0; j < n; j++ )
      {
        diag[j] = fmax ( diag[j], wa2[j] );
      }
    }
//
//  Beginning of the inner loop.
//
    for ( ; ; )
    {
//
//  Determine the direction P.
//
      dogleg ( n, r, lr, diag, qtf, delta, wa1, wa2, wa3 );
//
//  Store the direction P and X + P.  Calculate the norm of P.
//
      for ( j = 0; j < n; j++ )
      {
        wa1[j] = - wa1[j];
        wa2[j] = x[j] + wa1[j];
        wa3[j] = diag[j] * wa1[j];
      }
      pnorm = enorm ( n, wa3 );
//
//  On the first iteration, adjust the initial step bound.
//
      if ( iter == 1 )
      {
        delta = fmin ( delta, pnorm );
      }
//
//  Evaluate the function at X + P and calculate its norm.
//
      backward_euler_residual ( dydt, n, to, xo, t, wa2, wa4 );
      nfev = nfev + 1;
      fnorm1 = enorm ( n, wa4 );
//
//  Compute the scaled actual reduction.
//
      if ( fnorm1 < fnorm )
      {
        actred = 1.0 - ( fnorm1 / fnorm ) * ( fnorm1 / fnorm );
      }
      else
      {
        actred = - 1.0;
      }
//
//  Compute the scaled predicted reduction.
//
//  DO NOT ADJUST THIS LOOP, BECAUSE OF L.
//
      l = 1;
      for ( i = 1; i <= n; i++ )
      {
        sum = 0.0;
        for ( j = i; j <= n; j++ )
        {
          sum = sum + r[l-1] * wa1[j-1];
          l = l + 1;
        }
        wa3[i-1] = qtf[i-1] + sum;
      }
      temp = enorm ( n, wa3 );

      if ( temp < fnorm )
      {
        prered = 1.0 - ( temp / fnorm ) * ( temp / fnorm );
      }
      else
      {
        prered = 0.0;
      }
//
//  Compute the ratio of the actual to the predicted reduction.
//
      if ( 0.0 < prered )
      {
        ratio = actred / prered;
      }
      else
      {
        ratio = 0.0;
      }
//
//  Update the step bound.
//
      if ( ratio < p1 )
      {
        ncsuc = 0;
        ncfail = ncfail + 1;
        delta = p5 * delta;
      }
      else
      {
        ncfail = 0;
        ncsuc = ncsuc + 1;
        if ( p5 <= ratio || 1 < ncsuc )
        {
          delta = fmax ( delta, pnorm / p5 );
        }
        if ( fabs ( ratio - 1.0 ) <= p1 )
        {
          delta = pnorm / p5;
        }
      }
//
//  On successful iteration, update X, FVEC, and their norms.
//
      if ( p0001 <= ratio )
      {
        for ( j = 0; j < n; j++ )
        {
          x[j] = wa2[j];
          wa2[j] = diag[j] * x[j];
          fvec[j] = wa4[j];
        }
        xnorm = enorm ( n, wa2 );
        fnorm = fnorm1;
        iter = iter + 1;
      }
//
//  Determine the progress of the iteration.
//
      nslow1 = nslow1 + 1;
      if ( p001 <= actred )
      {
        nslow1 = 0;
      }
      if ( jeval )
      {
        nslow2 = nslow2 + 1;
      }
      if ( p1 <= actred )
      {
        nslow2 = 0;
      }
//
//  Test for convergence.
//
      if ( delta <= xtol * xnorm || fnorm == 0.0 )
      {
        info = 1;
        return info;
      }
//
//  Tests for termination and stringent tolerances.
//
      if ( maxfev <= nfev )
      {
        info = 2;
        return info;
      }
      if ( p1 * fmax ( p1 * delta, pnorm ) <= epsmch * xnorm )
      {
        info = 3;
        return info;
      }
      if ( nslow2 == 5 )
      {
        info = 4;
        return info;
      }
      if ( nslow1 == 10 )
      {
        info = 5;
        return info;
      }
//
//  Criterion for recalculating jacobian approximation by forward differences.
//
      if ( ncfail == 2 )
      {
        break;
      }
//
//  Calculate the rank one modification to the jacobian
//  and update QTF if necessary.
//
      for ( j = 0; j < n; j++ )
      {
        sum = 0.0;
        for ( i = 0; i < n; i++ )
        {
          sum = sum + fjac[i+j*ldfjac] * wa4[i];
        }
        wa2[j] = ( sum - wa3[j] ) / pnorm;
        wa1[j] = diag[j] * ( ( diag[j] * wa1[j] ) / pnorm );
        if ( p0001 <= ratio )
        {
          qtf[j] = sum;
        }
      }
//
//  Compute the QR factorization of the updated jacobian.
//
      r1updt ( n, n, r, lr, wa1, wa2, wa3 );
      r1mpyq ( n, n, fjac, ldfjac, wa2, wa3 );
      r1mpyq ( 1, n, qtf, 1, wa2, wa3 );

      jeval = false;
    }
//
//  End of the inner loop.
//
  }
//
//  End of the outer loop.
//
}
//****************************************************************************80

int hybrd_tr ( void dydt ( double t, double x[], double f[] ), int n, double to, 
  double xo[], double tn, double xn[],
  double fvec[], double xtol, int maxfev, int ml, int mu, double epsfcn,
  double diag[], int mode, double factor, int nfev,
  double fjac[], int ldfjac, double r[], int lr, double qtf[], double wa1[],
  double wa2[], double wa3[], double wa4[] )

//****************************************************************************80
//
//  Purpose:
//
//    hybrd_tr() finds a zero of a system of N nonlinear equations.
//
//  Discussion:
//
//    The code finds a zero of a system of
//    N nonlinear functions in N variables by a modification
//    of the Powell hybrid method. 
//
//    The jacobian is calculated by a forward-difference approximation.
//
//    The original code hybrd() was modified to deal with problems
//    involving a trapezoidal ODE solver residual.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    15 November 2023
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    This version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Input:
//
//    void dydt ( double t, double y[], double dy[] ):
//    the name of the user-supplied code which
//    evaluates the right hand side of the ODE.
//
//    int n: the number of functions and variables.
//
//    double to, xo[n]: the old time and solution.
// 
//    double tn, xn[n]: the new time and current solution estimate. 
//
//    double xtol: Termination occurs when the 
//    relative error between two consecutive iterates is at most XTOL.
//
//    int maxfev:  Termination occurs when the number of calls to FCN 
//    is at least MAXFEV by the end of an iteration.
//
//    int ml, mu: the number of subdiagonals and superdiagonals
//    within the band of the jacobian matrix.  If the jacobian is not 
//    banded, set ml and mu to at least n - 1.
//
//    double epsfcn: determines a suitable step length for the 
//    forward-difference approximation.  This approximation assumes that 
//    the relative errors in the functions are of the order of epsfcn. 
//    If epsfcn is less than the machine precision, it is assumed that 
//    the relative errors in the functions are of the order of the machine
//    precision.
//
//    double diag[n]: multiplicative scale factors for the 
//    variables.  Only needed as input if MODE = 2.
//
//    int mode: scaling option.
//    if mode = 1, the variables will be scaled internally. 
//    if mode = 2, the scaling is specified by the input diag.
//
//    double factor: determines the initial step bound. this bound is set 
//    to the product of factor and the euclidean norm of diag*x if nonzero, 
//    or else to factor itself. in most cases factor should lie in the
//    interval (.1,100.). 100. is a generally recommended value.
//
//    int ldfjac: not less than n, specifies the leading 
//    dimension of the array fjac.
//
//    int lr: the size of the R array, not less than (n*(n+1))/2.
//
//  Output:
//
//    double xn[n]: the final estimate of the solution.
//
//    double fvec[n]: the functions evaluated at the output value of X.
//
//    double diag[n]: multiplicative scale factors for the 
//    variables.  Set internally if MODE = 1.
//
//    int INFO: status flag.  A value of 1 indicates success. 
//    0: improper input parameters.
//    1: algorithm estimates that the relative error
//       between x and the solution is at most tol.
//    2: number of calls to the deriviative has reached or exceeded 200*(n+1).
//    3: tol is too small. no further improvement in
//       the approximate solution x is possible.
//    4: iteration is not making good progress.
//
//    int nfev: the number of calls to the derivative function.
//
//    double fjac[n*n]: the orthogonal matrix q produced by the qr 
//    factorization of the final approximate jacobian.
//
//    double r[lr]: the upper triangular matrix produced by the qr 
//    factorization of the final approximate jacobian, stored rowwise.
//
//    double qtf[n]: contains the vector (q transpose)*fvec.
//
//  Workspace:
//
//    double wa1[n], wa2[n], wa3[n], wa4[n].
//
{
  double actred;
  double delta;
  double epsmch;
  double fnorm;
  double fnorm1;
  int i;
  int info;
  int iter;
  int iwa[1];
  int j;
  bool jeval;
  int l;
  int msum;
  int ncfail;
  int ncsuc;
  int nslow1;
  int nslow2;
  const double p001 = 0.001;
  const double p0001 = 0.0001;
  const double p1 = 0.1;
  const double p5 = 0.5;
  double pnorm;
  double prered;
  double ratio;
  double sum;
  double temp;
  double xnorm;
//
//  Certain loops in this function were kept closer to their original FORTRAN77
//  format, to avoid confusing issues with the array index L.  These loops are
//  marked "DO NOT ADJUST", although they certainly could be adjusted (carefully)
//  once the initial translated code is tested.
//

//
//  EPSMCH is the machine precision.
//
  epsmch = std::numeric_limits<double>::epsilon();

  info = 0;
  nfev = 0;
//
//  Check the input parameters.
//
  if ( n <= 0 )
  {
    info = 0;
    return info;
  }
  if ( xtol < 0.0 )
  {
    info = 0;
    return info;
  }
  if ( maxfev <= 0 )
  {
    info = 0;
    return info;
  }
  if ( ml < 0 )
  {
    info = 0;
    return info;
  }
  if ( mu < 0 )
  {
    info = 0;
    return info;
  }
  if ( factor <= 0.0 )
  {
    info = 0;
    return info;
  }
  if ( ldfjac < n )
  {
    info = 0;
    return info;
  }
  if ( lr < ( n * ( n + 1 ) ) / 2 )
  {
    info = 0;
    return info;
  }
  if ( mode == 2 )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( diag[j] <= 0.0 )
      {
        info = 0;
        return info;
      }
    }
  }
//
//  Evaluate the function at the starting point and calculate its norm.
//
  trapezoidal_residual ( dydt, n, to, xo, tn, xn, fvec );
  nfev = 1;

  fnorm = enorm ( n, fvec );
//
//  Determine the number of calls to FCN needed to compute the jacobian matrix.
//
  if ( ml + mu + 1 < n )
  {
    msum = ml + mu + 1;
  }
  else
  {
    msum = n;
  }
//
//  Initialize iteration counter and monitors.
//
  iter = 1;
  ncsuc = 0;
  ncfail = 0;
  nslow1 = 0;
  nslow2 = 0;
//
//  Beginning of the outer loop.
//
  for ( ; ; )
  {
    jeval = true;
//
//  Calculate the jacobian matrix.
//
    fdjac_tr ( dydt, n, to, xo, tn, xn, fvec, fjac, ldfjac, ml, mu, epsfcn, 
      wa1, wa2 );

    nfev = nfev + msum;
//
//  Compute the QR factorization of the jacobian.
//
    qrfac ( n, n, fjac, ldfjac, false, iwa, 1, wa1, wa2 );
//
//  On the first iteration and if MODE is 1, scale according
//  to the norms of the columns of the initial jacobian.
//
    if ( iter == 1 )
    {
      if ( mode == 1 )
      {
        for ( j = 0; j < n; j++ )
        {
          if ( wa2[j] != 0.0 )
          {
            diag[j] = wa2[j];
          }
          else
          {
            diag[j] = 1.0;
          }
        }
      }
//
//  On the first iteration, calculate the norm of the scaled X
//  and initialize the step bound DELTA.
//
      for ( j = 0; j < n; j++ )
      {
        wa3[j] = diag[j] * xn[j];
      }
      xnorm = enorm ( n, wa3 );

      if ( xnorm == 0.0 )
      {
        delta = factor;
      }
      else
      {
        delta = factor * xnorm;
      }
    }
//
//  Form Q' * FVEC and store in QTF.
//
    for ( i = 0; i < n; i++ )
    {
      qtf[i] = fvec[i];
    }
    for ( j = 0; j < n; j++ )
    {
      if ( fjac[j+j*ldfjac] != 0.0 )
      {
        sum = 0.0;
        for ( i = j; i < n; i++ )
        {
          sum = sum + fjac[i+j*ldfjac] * qtf[i];
        }
        temp = - sum / fjac[j+j*ldfjac];
        for ( i = j; i < n; i++ )
        {
          qtf[i] = qtf[i] + fjac[i+j*ldfjac] * temp;
        }
      }
    }
//
//  Copy the triangular factor of the QR factorization into R.
//
//  DO NOT ADJUST THIS LOOP, BECAUSE OF L.
//
    for ( j = 1; j <= n; j++ )
    {
      l = j;
      for ( i = 1; i <= j - 1; i++ )
      {
        r[l-1] = fjac[(i-1)+(j-1)*ldfjac];
        l = l + n - i;
      }
      r[l-1] = wa1[j-1];
      if ( wa1[j-1] == 0.0 )
      {
        cerr << "  Matrix is singular.\n";
      }
    }
//
//  Accumulate the orthogonal factor in FJAC.
//
    qform ( n, n, fjac, ldfjac );
//
//  Rescale if necessary.
//
    if ( mode == 1 )
    {
      for ( j = 0; j < n; j++ )
      {
        diag[j] = fmax ( diag[j], wa2[j] );
      }
    }
//
//  Beginning of the inner loop.
//
    for ( ; ; )
    {
//
//  Determine the direction P.
//
      dogleg ( n, r, lr, diag, qtf, delta, wa1, wa2, wa3 );
//
//  Store the direction P and X + P.  Calculate the norm of P.
//
      for ( j = 0; j < n; j++ )
      {
        wa1[j] = - wa1[j];
        wa2[j] = xn[j] + wa1[j];
        wa3[j] = diag[j] * wa1[j];
      }
      pnorm = enorm ( n, wa3 );
//
//  On the first iteration, adjust the initial step bound.
//
      if ( iter == 1 )
      {
        delta = fmin ( delta, pnorm );
      }
//
//  Evaluate the function at X + P and calculate its norm.
//
      trapezoidal_residual ( dydt, n, to, xo, tn, wa2, wa4 );
      nfev = nfev + 1;
      fnorm1 = enorm ( n, wa4 );
//
//  Compute the scaled actual reduction.
//
      if ( fnorm1 < fnorm )
      {
        actred = 1.0 - ( fnorm1 / fnorm ) * ( fnorm1 / fnorm );
      }
      else
      {
        actred = - 1.0;
      }
//
//  Compute the scaled predicted reduction.
//
//  DO NOT ADJUST THIS LOOP, BECAUSE OF L.
//
      l = 1;
      for ( i = 1; i <= n; i++ )
      {
        sum = 0.0;
        for ( j = i; j <= n; j++ )
        {
          sum = sum + r[l-1] * wa1[j-1];
          l = l + 1;
        }
        wa3[i-1] = qtf[i-1] + sum;
      }
      temp = enorm ( n, wa3 );

      if ( temp < fnorm )
      {
        prered = 1.0 - ( temp / fnorm ) * ( temp / fnorm );
      }
      else
      {
        prered = 0.0;
      }
//
//  Compute the ratio of the actual to the predicted reduction.
//
      if ( 0.0 < prered )
      {
        ratio = actred / prered;
      }
      else
      {
        ratio = 0.0;
      }
//
//  Update the step bound.
//
      if ( ratio < p1 )
      {
        ncsuc = 0;
        ncfail = ncfail + 1;
        delta = p5 * delta;
      }
      else
      {
        ncfail = 0;
        ncsuc = ncsuc + 1;
        if ( p5 <= ratio || 1 < ncsuc )
        {
          delta = fmax ( delta, pnorm / p5 );
        }
        if ( fabs ( ratio - 1.0 ) <= p1 )
        {
          delta = pnorm / p5;
        }
      }
//
//  On successful iteration, update X, FVEC, and their norms.
//
      if ( p0001 <= ratio )
      {
        for ( j = 0; j < n; j++ )
        {
          xn[j] = wa2[j];
          wa2[j] = diag[j] * xn[j];
          fvec[j] = wa4[j];
        }
        xnorm = enorm ( n, wa2 );
        fnorm = fnorm1;
        iter = iter + 1;
      }
//
//  Determine the progress of the iteration.
//
      nslow1 = nslow1 + 1;
      if ( p001 <= actred )
      {
        nslow1 = 0;
      }
      if ( jeval )
      {
        nslow2 = nslow2 + 1;
      }
      if ( p1 <= actred )
      {
        nslow2 = 0;
      }
//
//  Test for convergence.
//
      if ( delta <= xtol * xnorm || fnorm == 0.0 )
      {
        info = 1;
        return info;
      }
//
//  Tests for termination and stringent tolerances.
//
      if ( maxfev <= nfev )
      {
        info = 2;
        return info;
      }
      if ( p1 * fmax ( p1 * delta, pnorm ) <= epsmch * xnorm )
      {
        info = 3;
        return info;
      }
      if ( nslow2 == 5 )
      {
        info = 4;
        return info;
      }
      if ( nslow1 == 10 )
      {
        info = 5;
        return info;
      }
//
//  Criterion for recalculating jacobian approximation by forward differences.
//
      if ( ncfail == 2 )
      {
        break;
      }
//
//  Calculate the rank one modification to the jacobian
//  and update QTF if necessary.
//
      for ( j = 0; j < n; j++ )
      {
        sum = 0.0;
        for ( i = 0; i < n; i++ )
        {
          sum = sum + fjac[i+j*ldfjac] * wa4[i];
        }
        wa2[j] = ( sum - wa3[j] ) / pnorm;
        wa1[j] = diag[j] * ( ( diag[j] * wa1[j] ) / pnorm );
        if ( p0001 <= ratio )
        {
          qtf[j] = sum;
        }
      }
//
//  Compute the QR factorization of the updated jacobian.
//
      r1updt ( n, n, r, lr, wa1, wa2, wa3 );
      r1mpyq ( n, n, fjac, ldfjac, wa2, wa3 );
      r1mpyq ( 1, n, qtf, 1, wa2, wa3 );

      jeval = false;
    }
//
//  End of the inner loop.
//
  }
//
//  End of the outer loop.
//
}
//****************************************************************************80

void qform ( int m, int n, double q[], int ldq )

//****************************************************************************80
//
//  Purpose:
//
//    qform() constructs the standard form of Q from its factored form.
//
//  Discussion:
//
//    This function proceeds from the computed QR factorization of
//    an M by N matrix A to accumulate the M by M orthogonal matrix
//    Q from its factored form.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    02 January 2018
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    This version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Input:
//
//    int M, the number of rows of A, and the order of Q.
//
//    int N, the number of columns of A.
//
//    double Q[LDQ*N].  The full lower trapezoid in
//    the first min(M,N) columns of Q contains the factored form.
//
//    int LDQ, the leading dimension of the array Q.
//
//  Output:
//
//    double Q[LDQ*N].  Q has been accumulated into a square matrix.
//
{
  int i;
  int j;
  int k;
  int minmn;
  double sum;
  double temp;
  double *wa;
//
//  Zero out the upper triangle of Q in the first min(M,N) columns.
//
  if ( m < n )
  {
    minmn = m;
  }
  else
  {
    minmn = n;
  }

  for ( j = 1; j < minmn; j++ )
  {
    for ( i = 0; i <= j - 1; i++ )
    {
      q[i+j*ldq] = 0.0;
    }
  }
//
//  Initialize remaining columns to those of the identity matrix.
//
  for ( j = n; j < m; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      q[i+j*ldq] = 0.0;
    }
    q[j+j*ldq] = 1.0;
  }
//
//  Accumulate Q from its factored form.
//
  wa = new double[m];

  for ( k = minmn - 1; 0 <= k; k-- )
  {
    for ( i = k; i < m; i++ )
    {
      wa[i] = q[i+k*ldq];
      q[i+k*ldq] = 0.0;
    }
    q[k+k*ldq] = 1.0;

    if ( wa[k] != 0.0 )
    {
      for ( j = k; j < m; j++ )
      {
        sum = 0.0;
        for ( i = k; i < m; i++ )
        {
          sum = sum + q[i+j*ldq] * wa[i];
        }
        temp = sum / wa[k];
        for ( i = k; i < m; i++ )
        {
          q[i+j*ldq] = q[i+j*ldq] - temp * wa[i];
        }
      }
    }
  }
//
//  Free memory.
//
  delete [] wa;

  return;
}
//****************************************************************************80

void qrfac ( int m, int n, double a[], int lda, bool pivot, int ipvt[],
  int lipvt, double rdiag[], double acnorm[] )

//****************************************************************************80
//
//  Purpose:
//
//    qrfac() computes the QR factorization of an M by N matrix.
//
//  Discussion:
//
//    This function uses Householder transformations with optional column
//    pivoting to compute a QR factorization of the M by N matrix A. 
//
//    That is, QRFAC determines an orthogonal
//    matrix Q, a permutation matrix P, and an upper trapezoidal
//    matrix R with diagonal elements of nonincreasing magnitude,
//    such that A*P = Q*R. 
//
//    The Householder transformation for
//    column k, k = 1,2,...,min(m,n), is of the form
//
//      i - (1/u(k))*u*u'
//
//    where U has zeros in the first K-1 positions. 
//
//    The form of this transformation and the method of pivoting first
//    appeared in the corresponding LINPACK function.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    02 January 2017
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    This version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Parameters:
//
//    Input, int M, the number of rows of A.
//
//    Input, int N, the number of columns of A.
//
//    Input/output, double A[M*N].  On input, the matrix for which the QR 
//    factorization is to be computed.  On output, the strict upper trapezoidal 
//    part contains the strict upper trapezoidal part of the R factor, and 
//    the lower trapezoidal part contains a factored form of Q, the non-trivial
//    elements of the U vectors described above.
//
//    Input, int LDA, a positive value not less than M which specifies the 
//    leading dimension of the array A.
//
//    Input, bool PIVOT.  If true, then column pivoting is enforced.
//
//    Output, int IPVT[LIPVT].  If PIVOT is true, then on output IPVT 
//    defines the permutation matrix P such that A*P = Q*R.  Column J of P
//    is column IPVT[J] of the identity matrix.
//
//    Input, int LIPVT.  If PIVOT is false, then LIPVT may be as small 
//    as 1.  If PIVOT is true, then LIPVT must be at least N.
//
//       rdiag is an output array of length n which contains the
//         diagonal elements of r.
//
//       acnorm is an output array of length n which contains the
//         norms of the corresponding columns of the input matrix a.
//         if this information is not needed, then acnorm can coincide
//         with rdiag.
//
{
  double ajnorm;
  double epsmch;
  int i;
  int j;
  int k;
  int kmax;
  int minmn;
  const double p05 = 0.05;
  double sum;
  double temp;
  double *wa;
//
//  EPSMCH is the machine precision.
//
  epsmch = std::numeric_limits<double>::epsilon();
//
//  Compute the initial column norms and initialize several arrays.
//
  wa = new double[n];

  for ( j = 0; j < n; j++ )
  {
    acnorm[j] = enorm ( m, a+j*lda );
    rdiag[j] = acnorm[j];
    wa[j] = rdiag[j];
    if ( pivot )
    {
      ipvt[j] = j;
    }
  }
//
//  Reduce A to R with Householder transformations.
//
  if ( m < n )
  {
    minmn = m;
  }
  else
  {
    minmn = n;
  }

  for ( j = 0; j < minmn; j++ )
  {
    if ( pivot )
    {
//
//  Bring the column of largest norm into the pivot position.
//
      kmax = j;
      for ( k = j; k < n; k++ )
      {
        if ( rdiag[kmax] < rdiag[k] )
        {
          kmax = k;
        }
      }
      if ( kmax != j )
      {
        for ( i = 0; i < m; i++ )
        {
          temp          = a[i+j*lda];
          a[i+j*lda]    = a[i+kmax*lda];
          a[i+kmax*lda] = temp;
        }
        rdiag[kmax] = rdiag[j];
        wa[kmax]    = wa[j];
        k          = ipvt[j];
        ipvt[j]    = ipvt[kmax];
        ipvt[kmax] = k;
      }
    }
//
//  Compute the Householder transformation to reduce the
//  J-th column of A to a multiple of the J-th unit vector.
//
    ajnorm = enorm ( m - j, a+j+j*lda );

    if ( ajnorm != 0.0 )
    {
      if ( a[j+j*lda] < 0.0 )
      {
        ajnorm = - ajnorm;
      }
      for ( i = j; i < m; i++ )
      {
        a[i+j*lda] = a[i+j*lda] / ajnorm;
      }
      a[j+j*lda] = a[j+j*lda] + 1.0;
//
//  Apply the transformation to the remaining columns and update the norms.
//
      for ( k = j + 1; k < n; k++ )
      {
        sum = 0.0;
        for ( i = j; i < m; i++ )
        {
          sum = sum + a[i+j*lda] * a[i+k*lda];
        }
        temp = sum / a[j+j*lda];
        for ( i = j; i < m; i++ )
        {
          a[i+k*lda] = a[i+k*lda] - temp * a[i+j*lda];
        }
        if ( pivot && rdiag[k] != 0.0 )
        {
          temp = a[j+k*lda] / rdiag[k];
          rdiag[k] = rdiag[k] * sqrt ( fmax ( 0.0, 1.0 - temp * temp ) );
          if ( p05 * ( rdiag[k] / wa[k] ) * ( rdiag[k] / wa[k] ) <= epsmch )
          {
            rdiag[k] = enorm ( m - 1 - j, a+(j+1)+k*lda );
            wa[k] = rdiag[k];
          }
        }
      }
    }
    rdiag[j] = - ajnorm;
  }
//
//  Free memory.
//
  delete [] wa;

  return;
}
//****************************************************************************80

void r1mpyq ( int m, int n, double a[], int lda, double v[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    r1mpyq() multiplies an M by N matrix A by the Q factor.
//
//  Discussion:
//
//    Given an M by N matrix A, this function computes a*q where
//    q is the product of 2*(n - 1) transformations
//
//      gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
//
//    and gv(i), gw(i) are givens rotations in the (i,n) plane which
//    eliminate elements in the i-th and n-th planes, respectively.
//
//    Q itself is not given, rather the information to recover the
//    GV and GW rotations is supplied.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    05 April 2010
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    This version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Input:
//
//    int M, the number of rows of A.
//
//    int N, the number of columns of A.
//
//    double A[M*N].  The matrix to be postmultiplied  by Q.
//
//    int LDA, a positive value not less than M
//    which specifies the leading dimension of the array A.
//
//    double V[N].  V(I) must contain the information necessary to 
//    recover the givens rotation GV(I) described above.
//
//    double W[N], contains the information necessary to recover the
//    Givens rotation gw(i) described above.
//
//  Output:
//
//    double A[M*N].  The value of A*Q.
//
{
  double c;
  int i;
  int j;
  double s;
  double temp;
//
//  Apply the first set of Givens rotations to A.
//
  for ( j = n - 2; 0 <= j; j-- )
  {
    if ( 1.0 < fabs ( v[j] ) )
    {
      c = 1.0 / v[j];
      s = sqrt ( 1.0 - c * c );
    }
    else
    {
      s = v[j];
      c = sqrt ( 1.0 - s * s );
    }
    for ( i = 0; i < m; i++ )
    {
      temp           = c * a[i+j*lda] - s * a[i+(n-1)*lda];
      a[i+(n-1)*lda] = s * a[i+j*lda] + c * a[i+(n-1)*lda];
      a[i+j*lda]     = temp;
    }
  }
//
//  Apply the second set of Givens rotations to A.
//
  for ( j = 0; j < n - 1; j++ )
  {
    if ( 1.0 < fabs ( w[j] ) )
    {
      c = 1.0 / w[j];
      s = sqrt ( 1.0 - c * c );
    }
    else
    {
      s = w[j];
      c = sqrt ( 1.0 - s * s );
    }
    for ( i = 0; i < m; i++ )
    {
      temp           =   c * a[i+j*lda] + s * a[i+(n-1)*lda];
      a[i+(n-1)*lda] = - s * a[i+j*lda] + c * a[i+(n-1)*lda];
      a[i+j*lda]     = temp;
    }
  }

  return;
}
//****************************************************************************80

bool r1updt ( int m, int n, double s[], int ls, double u[], double v[],
  double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    r1updt() updates the Q factor after a rank one update of the matrix.
//
//  Discussion:
//
//    Given an M by N lower trapezoidal matrix S, an M-vector U,
//    and an N-vector V, the problem is to determine an
//    orthogonal matrix Q such that
//
//      (S + U*V') * Q
//
//    is again lower trapezoidal.
//
//    This function determines q as the product of 2*(n - 1) transformations
//
//      gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
//
//    where gv(i), gw(i) are givens rotations in the (i,n) plane
//    which eliminate elements in the i-th and n-th planes,
//    respectively. 
//
//    Q itself is not accumulated, rather the
//    information to recover the gv, gw rotations is returned.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    08 April 2010
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    This version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Parameters:
//
//    Input, int M, the number of rows of S.
//
//    Input, int N, the number of columns of S.  N must not exceed M.
//
//    Input/output, double S[LS].  On input, the lower trapezoidal matrix S 
//    stored by columns.  On output, the lower trapezoidal matrix produced 
//    as described above.
//
//    Input, int LS, a positive value not less than (N*(2*M-N+1))/2.
//
//    Input, double U[M], contains the vector U.
//
//    Input/output, double V[N].  On input, the vector v.  On output,
//    information necessary to recover the givens rotation gv(i) described above.
//
//    double W[M]: information necessary to recover the givens rotation gw(i).
//
//    bool SING: true if any diagonal elements of S is zero.
//
{
  double cotan;
  double cs;
  double giant;
  int i;
  int j;
  int jj;
  int l;
  int nm1;
  const double p25 = 0.25;
  const double p5 = 0.5;
  double sn;
  bool sing;
  double tan;
  double tau;
  double temp;
//
//  Because of the computation of the pointer JJ, this function was
//  converted from FORTRAN77 to C++ in a conservative way.  All computations
//  are the same, and only array indexing is adjusted.
//
//  GIANT is the largest magnitude.
//
  giant = std::numeric_limits<double>::max();
//
//  Initialize the diagonal element pointer.
//
  jj = ( n * ( 2 * m - n + 1 ) ) / 2 - ( m - n );
//
//  Move the nontrivial part of the last column of S into W.
//
  l = jj;
  for ( i = n; i <= m; i++ )
  {
    w[i-1] = s[l-1];
    l = l + 1;
  }
//
//  Rotate the vector V into a multiple of the N-th unit vector
//  in such a way that a spike is introduced into W.
//
  nm1 = n - 1;

  for ( j = n - 1; 1 <= j; j-- )
  {
    jj = jj - ( m - j + 1 );
    w[j-1] = 0.0;

    if ( v[j-1] != 0.0 )
    {
//
//  Determine a Givens rotation which eliminates the J-th element of V.
//
      if ( fabs ( v[n-1] ) < fabs ( v[j-1] ) )
      {
        cotan = v[n-1] / v[j-1];
        sn = p5 / sqrt ( p25 + p25 * cotan * cotan );
        cs = sn * cotan;
        tau = 1.0;
        if ( 1.0 < fabs ( cs ) * giant )
        {
          tau = 1.0 / cs;
        }
      }
      else
      {
        tan = v[j-1] / v[n-1];
        cs = p5 / sqrt ( p25 + p25 * tan * tan );
        sn = cs * tan;
        tau = sn;
      }
//
//  Apply the transformation to V and store the information
//  necessary to recover the Givens rotation.
//
      v[n-1] = sn * v[j-1] + cs * v[n-1];
      v[j-1] = tau;
//
//  Apply the transformation to S and extend the spike in W.
//
      l = jj;
      for ( i = j; i <= m; i++ )
      {
        temp   = cs * s[l-1] - sn * w[i-1];
        w[i-1] = sn * s[l-1] + cs * w[i-1];
        s[l-1] = temp;
        l = l + 1;
      }
    }
  }
//
//  Add the spike from the rank 1 update to W.
//
  for ( i = 1; i <= m; i++ )
  {
     w[i-1] = w[i-1] + v[n-1] * u[i-1];
  }
//
//  Eliminate the spike.
//
  sing = false;

  for ( j = 1; j <= nm1; j++ )
  {
//
//  Determine a Givens rotation which eliminates the
//  J-th element of the spike.
//
    if ( w[j-1] != 0.0 )
    {

      if ( fabs ( s[jj-1] ) < fabs ( w[j-1] ) )
      {
        cotan = s[jj-1] / w[j-1];
        sn = p5 / sqrt ( p25 + p25 * cotan * cotan );
        cs = sn * cotan;
        tau = 1.0;
        if ( 1.0 < fabs ( cs ) * giant )
        {
          tau = 1.0 / cs;
        }
      }
      else
      {
        tan = w[j-1] / s[jj-1];
        cs = p5 / sqrt ( p25 + p25 * tan * tan );
        sn = cs * tan;
        tau = sn;
      }
//
//  Apply the transformation to s and reduce the spike in w.
//
      l = jj;

      for ( i = j; i <= m; i++ )
      {
        temp   =   cs * s[l-1] + sn * w[i-1];
        w[i-1] = - sn * s[l-1] + cs * w[i-1];
        s[l-1] = temp;
        l = l + 1;
      }
//
//  Store the information necessary to recover the givens rotation.
//
      w[j-1] = tau;
    }
//
//  Test for zero diagonal elements in the output s.
//
    if ( s[jj-1] == 0.0 )
    {
      sing = true;
    }
    jj = jj + ( m - j + 1 );
  }
//
//  Move W back into the last column of the output S.
//
  l = jj;
  for ( i = n; i <= m; i++ )
  {
    s[l-1] = w[i-1];
    l = l + 1;
  }
  if ( s[jj-1] == 0.0 )
  {
    sing = true;
  }
  return sing;
}
//****************************************************************************80

void trapezoidal_residual ( void dydt ( double t, double y[], double dy[] ), 
  int n, double to, double yo[], double tn, double yn[], double ft[] )

//****************************************************************************80
//
//  Purpose:
//
//    trapezoidal_residual() evaluates the trapezoidal ODE residual.
//
//  Discussion:
//
//    Let to and tn be two times, with yo and yn the associated ODE
//    solution values there.  If yn satisfies the trapezoidal ODE condition,
//    then
//
//      0.5 * ( dydt(to,yo) + dydt(tn,yn) ) = ( yn - yo ) / ( tn - to )
//
//    This can be rewritten as
//
//      residual = yn - yo - ( tn - to ) * 0.5 * ( dydt(to,yo) + dydt(tn,yn) )
//
//    Given the other information, a nonlinear equation solver can be used
//    to estimate the value yn that makes the residual zero.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    16 November 2023
//
//  Author:
//
//    John Burkardt
//
//  Input:
//
//    void dydt( double t, double y[], double dy[] ): 
//    evaluates the right hand side of the ODE.
//
//    int n: the vector size.
//
//    double to, yo[n]: the old time and solution.
//
//    double tn, yn[n]: the new time and tentative solution.
//
//  Output:
//
//    double ft[n]: the residual.
//
{
  int i;
  double *dydtn;
  double *dydto;
  
  dydtn = new double[n];
  dydto = new double[n];

  dydt ( tn, yn, dydtn );
  dydt ( to, yo, dydto );

  for ( i = 0; i < n; i++ )
  {
    ft[i] = yn[i] - yo[i] - ( tn - to ) * 0.5 * ( dydtn[i] + dydto[i] );
  }

  delete [] dydtn;
  delete [] dydto;

  return;
}
