void backward_euler_residual ( void dydt ( double t, double y[], double dy[] ), 
  int n, double to, double yo[], double tm, double ym[], double fm[] );
void bdf2_residual ( void dydt ( double t, double y[], double dy[] ), 
  int n, double t1, double y1[], double t2, double y2[], double t3, double y3[], 
  double fm[] );
void dogleg ( int n, double r[], int lr, double diag[], double qtb[],
  double delta, double x[], double wa1[], double wa2[] );
double enorm ( int n, double x[] );
void fdjac1 ( void fcn ( int n, double x[], double f[] ), 
  int n, double x[], double fvec[], double fjac[], int ldfjac,
  int ml, int mu, double epsfcn, double wa1[], double wa2[] );
void fdjac_bdf2 ( void dydt ( double t, double x[], double f[] ),
  int n, double t1, double x1[], double t2, double x2[], double t3, double x3[], 
  double fvec[], double fjac[], int ldfjac, int ml, int mu, double epsfcn, 
  double wa1[], double wa2[] );
void fdjac_be ( void dydt ( double t, double x[], double f[] ),
  int n, double to, double xo[], double t, double x[], double fvec[], 
  double fjac[], int ldfjac, int ml, int mu, double epsfcn, double wa1[], 
  double wa2[] );
void fdjac_tr ( void dydt ( double t, double x[], double f[] ),
  int n, double to, double xo[], double tn, double xn[], double fvec[], 
  double fjac[], int ldfjac, int ml, int mu, double epsfcn, double wa1[], 
  double wa2[] );
int fsolve ( void fcn ( int n, double x[], double fvec[] ), int n, 
  double x[], double fvec[], double tol, double wa[], int lwa );
int fsolve_bdf2 ( void dydt ( double t, double x[], double f[] ), int n, 
  double t1, double x1[], double t2, double x2[], double t3, double x3[], 
  double fvec[], double tol );
int fsolve_be ( void dydt ( double t, double x[], double f[] ), int n, 
  double to, double xo[], double t, double x[], double fvec[], double tol );
int fsolve_tr ( void dydt ( double t, double x[], double f[] ), int n, 
  double to, double xo[], double tn, double xn[], double fvec[], double tol );
int hybrd ( void fcn ( int n, double x[], double fvec[] ), 
  int n, double x[], double fvec[], double xtol, int maxfev, int ml, 
  int mu, double epsfcn, double diag[], int mode, double factor, 
  int nfev, double fjac[], int ldfjac, double r[], int lr, double qtf[], 
  double wa1[], double wa2[], double wa3[], double wa4[] );
int hybrd_bdf2 ( void dydt ( double t, double x[], double f[] ), int n, 
  double t1, double x1[], double t2, double x2[], double t3, double x3[],
  double fvec[], double xtol, int maxfev, int ml, int mu, double epsfcn,
  double diag[], int mode, double factor, int nfev,
  double fjac[], int ldfjac, double r[], int lr, double qtf[], double wa1[],
  double wa2[], double wa3[], double wa4[] );
int hybrd_be ( void dydt ( double t, double x[], double f[] ), int n, double to, 
  double xo[], double t, double x[],
  double fvec[], double xtol, int maxfev, int ml, int mu, double epsfcn,
  double diag[], int mode, double factor, int nfev,
  double fjac[], int ldfjac, double r[], int lr, double qtf[], double wa1[],
  double wa2[], double wa3[], double wa4[] );
int hybrd_tr ( void dydt ( double t, double x[], double f[] ), int n, double to, 
  double xo[], double tn, double xn[],
  double fvec[], double xtol, int maxfev, int ml, int mu, double epsfcn,
  double diag[], int mode, double factor, int nfev,
  double fjac[], int ldfjac, double r[], int lr, double qtf[], double wa1[],
  double wa2[], double wa3[], double wa4[] );
void qform ( int m, int n, double q[], int ldq );
void qrfac ( int m, int n, double a[], int lda, bool pivot, int ipvt[],
  int lipvt, double rdiag[], double acnorm[] );
void r1mpyq ( int m, int n, double a[], int lda, double v[], double w[] );
bool r1updt ( int m, int n, double s[], int ls, double u[], double v[], double w[] );
void trapezoidal_residual ( void dydt ( double t, double y[], double dy[] ), 
  int n, double to, double yo[], double tn, double yn[], double ft[] );
