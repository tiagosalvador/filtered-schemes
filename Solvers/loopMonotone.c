#include "mex.h"
#include "math.h"

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

/* function to convert indexes in Matlab to indexes in c */
mwSize IND(mwSize i, mwSize j, mwSize n){
    return i+n*j;
}

/* function to check if x is a member of the matrix S */
int isMember(double x, double y, double *S, mwSize n){
   int r=0;
   mwSize i;
   for (i=0;i<n;i++){
       if (S[IND(i,0,n)]==x && S[IND(i,1,n)]==y){r++;}
   }
   return r;
}

void update(mwSize i, mwSize j, double *u, mwSize n, double *f, double *x, double *y){
  double a,b,dx;
  mwSize i1,i2,j1,j2;
  dx = x[1]-x[0];
  if (i-1 < 0)
      a = u[IND(i+1,j,n)];
  else if (i+1 > n-1)
      a = u[IND(i-1,j,n)];
  else
      a = MIN(u[IND(i+1,j,n)],u[IND(i-1,j,n)]);
  if (j-1 < 0)
      b = u[IND(i,j+1,n)];
  else if (j+1 > n-1)
      b = u[IND(i,j-1,n)];
  else
      b = MIN(u[IND(i,j+1,n)],u[IND(i,j-1,n)]);
  if (MAX(a-b,b-a) > f[IND(i,j,n)]*dx)
      u[IND(i,j,n)] = MIN(a,b) + f[IND(i,j,n)]*dx;
  else
      u[IND(i,j,n)] = (a+b+sqrt(2*pow(f[IND(i,j,n)]*dx,2)-pow(a-b,2)))/2;
}

void loopSweeping(double *u, double *f, mwSize n, double *dp, mwSize ndp, double *x, double *y){
    
    mwSize i,j;
    
    /* loop left, right, up, down */
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            if (dp[IND(i,j,ndp)] == 0){
                update(i,j,u,n,f,x,y);
            }
        }
    }
    /* loop right, left, up, down */
    for (i=n-1; i>-1; i--) {
        for (j=0; j<n; j++) {
            if (dp[IND(i,j,ndp)] == 0){
                update(i,j,u,n,f,x,y);
            }
        }
    }
    /* loop left, right, down, up */
    for (i=0; i<n; i++) {
        for (j=n-1; j>-1; j--) {
            if (dp[IND(i,j,ndp)] == 0){
                update(i,j,u,n,f,x,y);
            }
        }
    }
    /* loop right, left, down, up */
    for (i=n-1; i>-1; i--) {
        for (j=n-1; j>-1; j--) {
            if (dp[IND(i,j,ndp)] == 0){
                update(i,j,u,n,f,x,y);
            }
        }
    }
    
    /* loop up, down, left, right */
    for (j=0; j<n; j++) {
        for (i=0; i<n; i++) {
            if (dp[IND(i,j,ndp)] == 0){
                update(i,j,u,n,f,x,y);
            }
        }
    }
    /* loop up, down, right, left */
    for (j=0; j<n; j++) {
        for (i=n-1; i>-1; i--) {
            if (dp[IND(i,j,ndp)] == 0){
                update(i,j,u,n,f,x,y);
            }
        }
    }
    /* loop down, up, left, right */
    for (j=n-1; j>-1; j--) {
        for (i=0; i<n; i++) {
            if (dp[IND(i,j,ndp)] == 0){
                update(i,j,u,n,f,x,y);
            }
        }
    }
    /* loop down, up, right, left */
    for (j=n-1; j>-1; j--) {
        for (i=n-1; i>-1; i--) {
            if (dp[IND(i,j,ndp)] == 0){
                update(i,j,u,n,f,x,y);
            }
        }
    }
}



void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *u, *f, *dp, *x, *y, *unew;
  mwSize n,ndp,i,j,count;
  
  
  /* obtain the pointers to the input data */
  u = mxGetPr(prhs[0]);
  f = mxGetPr(prhs[1]);
  dp = mxGetPr(prhs[2]);
  x = mxGetPr(prhs[3]);
  y = mxGetPr(prhs[4]);
  
  n = mxGetM(prhs[0]);
  ndp = mxGetM(prhs[2]);
  
  /*  call the C subroutine */
  loopSweeping(u,f,n,dp,ndp,x,y);
  
  plhs[0] = mxCreateDoubleMatrix(n,n, mxREAL);
  
  unew = mxGetPr(plhs[0]);
  
  /* copy the updated variable u and argmin to the output variable */
  count = 0;
  for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
          unew[count] = u[count];
          count++;
      }
  }
  
}