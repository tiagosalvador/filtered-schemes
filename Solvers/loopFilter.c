#include "mex.h"
#include "math.h"

// #if !defined(MAX)
// #define	MAX(A, B)	((A) > (B) ? (A) : (B))
// #endif

// #if !defined(MIN)
// #define	MIN(A, B)	((A) < (B) ? (A) : (B))
// #endif

/* function to convert indexes in Matlab to indexes in c */
mwSize IND(mwSize i, mwSize j, mwSize n){
    return i+n*j;
}

double MIN(double a, double b){
    double r;
    if (a <= b)
        r = a;
    else
        r = b;
    return r;
}

double MAX(double a, double b){
    double r;
    if (a >= b)
        r = a;
    else
        r = b;
    return r;
}

/* function to check if x is a member of the matrix S */
mwSize isMember(double x, double y, double *S, mwSize n){
   mwSize r=0;
   mwSize i;
   for (i=0;i<n;i++){
       if (S[IND(i,0,n)]==x && S[IND(i,1,n)]==y){r = r+1;}
   }
   return r;
}

double absolutevalue(double x){
    return MAX(x,-x);
}

double solve(double a, double b, double c){
    double r;
    if (absolutevalue(a-b) > c)
        r = MIN(a,b)+c;
    else
        r = (a+b+sqrt(2*pow(c,2)-pow(a-b,2)))/2;
    return r;
}

double Fmonotone(mwSize i, mwSize j, double *u, mwSize n, double *f, double *x, double *y){
    double uij,uL,uR,uT,uB,dx,dy,r;
    mwSize i1,i2,j1,j2;
    dx = x[1]-x[0];
    dy = y[1]-y[0];
    i1 = MAX(i-1,0);
    i2 = MIN(n-1,i+1);
    j1 = MAX(j-1,0);
    j2 = MIN(n-1,j+1);
    uij = u[IND(i,j,n)];
    uT = u[IND(i1,j,n)];
    uB = u[IND(i2,j,n)];
    uL = u[IND(i,j1,n)];
    uR = u[IND(i,j2,n)];
    r = sqrt(pow(MAX(MAX(uij-uL,uij-uR),0)/dx,2) + pow(MAX(MAX(uij-uT,uij-uB),0)/dy,2))-f[IND(i,j,n)];
    return r;
}

void Fmonotonefast(mwSize i, mwSize j, double *u, mwSize n, double *f, double *x, double *y){
  double a,b,c,dx;
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
  c = f[IND(i,j,n)]*dx;
  u[IND(i,j,n)] = MIN(solve(a,b,c),u[IND(i,j,n)]);
}



double Faccurate2nd(mwSize i, mwSize j, double *u, mwSize n, double *f, double *x, double *y){
    double rij,dx,dy,ux,uy,uij,uL,uR,uT,uB;
    mwSize i1,i2,j1,j2;
    dx = x[1]-x[0];
    dy = y[1]-y[0];
    i1 = MAX(i-1,0);
    i2 = MIN(n-1,i+1);
    j1 = MAX(j-1,0);
    j2 = MIN(n-1,j+1);
    uij = u[IND(i,j,n)];
    uT = u[IND(i1,j,n)];
    uB = u[IND(i2,j,n)];
    uL = u[IND(i,j1,n)];
    uR = u[IND(i,j2,n)];
    if (j>0 && j<n-1)
        ux = (uL-uR)/(2*dx);
    else
        ux = MAX(MAX(uij-uL,uij-uR),0)/dx;
    if (i>0 && i<n-1)
        uy = (uT-uB)/(2*dy);
    else
        uy = MAX(MAX(uij-uT,uij-uB),0)/dy;
    rij = sqrt(pow(ux,2)+pow(uy,2))-f[IND(i,j,n)];
    return rij;
}

double Faccurate2ndupwind(mwSize i, mwSize j, double *u, mwSize n, double *f, double *x, double *y){
    double rij,dx,dy,ux,uy;
    dx = x[1]-x[0];
    dy = y[1]-y[0];
    if (j < 2)
        ux = MAX(3*u[IND(i,j,n)]-4*u[IND(i,j+1,n)]+u[IND(i,j+2,n)],0)/(2*dx);
    else if (j > n-3)
        ux = MAX(3*u[IND(i,j,n)]-4*u[IND(i,j-1,n)]+u[IND(i,j-2,n)],0)/(2*dx);
    else
        ux = MAX(MAX(3*u[IND(i,j,n)]-4*u[IND(i,j+1,n)]+u[IND(i,j+2,n)],3*u[IND(i,j,n)]-4*u[IND(i,j-1,n)]+u[IND(i,j-2,n)]),0)/(2*dx);
    if (i < 2)
        uy = MAX(3*u[IND(i,j,n)]-4*u[IND(i+1,j,n)]+u[IND(i+2,j,n)],0)/(2*dx);
    else if (i > n-3)
        uy = MAX(3*u[IND(i,j,n)]-4*u[IND(i-1,j,n)]+u[IND(i-2,j,n)],0)/(2*dx);
    else
        uy = MAX(MAX(3*u[IND(i,j,n)]-4*u[IND(i+1,j,n)]+u[IND(i+2,j,n)],3*u[IND(i,j,n)]-4*u[IND(i-1,j,n)]+u[IND(i-2,j,n)]),0)/(2*dx);
    rij = sqrt(pow(ux,2)+pow(uy,2))-f[IND(i,j,n)];
    return rij;
}
    
void Faccurate2ndupwindfast(mwSize i, mwSize j, double *u, mwSize n, double *f, double *x, double *y){
    double dx,a,b,c;
    dx = x[1]-x[0];
    if (i-2 < 0)
        a = (4*u[IND(i+1,j,n)]-u[IND(i+2,j,n)])/3;
    else if (i+2 > n-1)
        a = (4*u[IND(i-1,j,n)]-u[IND(i-2,j,n)])/3;
    else
        a = MIN(4*u[IND(i+1,j,n)]-u[IND(i+2,j,n)],4*u[IND(i-1,j,n)]-u[IND(i-2,j,n)])/3;
    if (j-2 < 0)
        b = (4*u[IND(i,j+1,n)]-u[IND(i,j+2,n)])/3;
    else if (j+2 > n-1)
        b = (4*u[IND(i,j-1,n)]-u[IND(i,j-2,n)])/3;
    else
        b = MIN(4*u[IND(i,j+1,n)]-u[IND(i,j+2,n)],4*u[IND(i,j-1,n)]-u[IND(i,j-2,n)])/3;
    c = f[IND(i,j,n)]*2*dx/3;
    u[IND(i,j,n)] = solve(a,b,c);
}

double Faccurate3rdupwind(mwSize i, mwSize j, double *u, mwSize n, double *f, double *x, double *y){
    double rij,dx,dy,ux,uy;
    dx = x[1]-x[0];
    dy = y[1]-y[0];
    if (j-3 < 0)
        ux = MAX(11*u[IND(i,j,n)]-18*u[IND(i,j+1,n)]+9*u[IND(i,j+2,n)]-2*u[IND(i,j+3,n)],0)/(6*dx);
    else if (j+3 > n-1)
        ux = MAX(11*u[IND(i,j,n)]-18*u[IND(i,j-1,n)]+9*u[IND(i,j-2,n)]-2*u[IND(i,j-3,n)],0)/(6*dx);
    else
        ux = MAX(MAX(11*u[IND(i,j,n)]-18*u[IND(i,j+1,n)]+9*u[IND(i,j+2,n)]-2*u[IND(i,j+3,n)],11*u[IND(i,j,n)]-18*u[IND(i,j-1,n)]+9*u[IND(i,j-2,n)]-2*u[IND(i,j-3,n)]),0)/(6*dx);
    if (i-3 < 0)
        uy = MAX(11*u[IND(i,j,n)]-18*u[IND(i+1,j,n)]+9*u[IND(i+2,j,n)]-2*u[IND(i+3,j,n)],0)/(6*dx);
    else if (i+3 > n-1)
        uy = MAX(11*u[IND(i,j,n)]-18*u[IND(i-1,j,n)]+9*u[IND(i-2,j,n)]-2*u[IND(i-3,j,n)],0)/(6*dx);
    else
        uy = MAX(MAX(11*u[IND(i,j,n)]-18*u[IND(i+1,j,n)]+9*u[IND(i+2,j,n)]-2*u[IND(i+3,j,n)],11*u[IND(i,j,n)]-18*u[IND(i-1,j,n)]+9*u[IND(i-2,j,n)]-2*u[IND(i-3,j,n)]),0)/(6*dx);
    rij = sqrt(pow(ux,2)+pow(uy,2))-f[IND(i,j,n)];
    return rij;
}

void Faccurate3rdupwindfast(mwSize i, mwSize j, double *u, mwSize n, double *f, double *x, double *y){
    double dx,a,b,c;
    dx = x[1]-x[0];
    if (i-3 < 0)
        a = (18*u[IND(i+1,j,n)]-9*u[IND(i+2,j,n)]+2*u[IND(i+3,j,n)])/11;
    else if (i+3 > n-1)
        a = (18*u[IND(i-1,j,n)]-9*u[IND(i-2,j,n)]+2*u[IND(i-3,j,n)])/11;
    else
        a = MIN(18*u[IND(i+1,j,n)]-9*u[IND(i+2,j,n)]+2*u[IND(i+3,j,n)],18*u[IND(i-1,j,n)]-9*u[IND(i-2,j,n)]+2*u[IND(i-3,j,n)])/11;
    if (j-3 < 0)
        b = (18*u[IND(i,j+1,n)]-9*u[IND(i,j+2,n)]+2*u[IND(i,j+3,n)])/11;
    else if (j+3 > n-1)
        b = (18*u[IND(i,j-1,n)]-9*u[IND(i,j-2,n)]+2*u[IND(i,j-3,n)])/11;
    else
        b = MIN(18*u[IND(i,j+1,n)]-9*u[IND(i,j+2,n)]+2*u[IND(i,j+3,n)],18*u[IND(i,j-1,n)]-9*u[IND(i,j-2,n)]+2*u[IND(i,j-3,n)])/11;
    c = f[IND(i,j,n)]*6*dx/11;
    u[IND(i,j,n)] = solve(a,b,c);
}

double Faccurate4thupwind(mwSize i, mwSize j, double *u, mwSize n, double *f, double *x, double *y){
    double rij,dx,dy,ux,uy;
    dx = x[1]-x[0];
    dy = y[1]-y[0];
    if (j-4 < 0)
        ux = MAX(25*u[IND(i,j,n)]-48*u[IND(i,j+1,n)]+36*u[IND(i,j+2,n)]-16*u[IND(i,j+3,n)]+3*u[IND(i,j+4,n)],0)/(12*dx);
    else if (j+4 > n-1)
        ux = MAX(25*u[IND(i,j,n)]-48*u[IND(i,j-1,n)]+36*u[IND(i,j-2,n)]-16*u[IND(i,j-3,n)]+3*u[IND(i,j-4,n)],0)/(12*dx);
    else
        ux = MAX(MAX(25*u[IND(i,j,n)]-48*u[IND(i,j+1,n)]+36*u[IND(i,j+2,n)]-16*u[IND(i,j+3,n)]+3*u[IND(i,j+4,n)],25*u[IND(i,j,n)]-48*u[IND(i,j-1,n)]+36*u[IND(i,j-2,n)]-16*u[IND(i,j-3,n)]+3*u[IND(i,j-4,n)]),0)/(12*dx);
    if (i-4 < 0)
        uy = MAX(25*u[IND(i,j,n)]-48*u[IND(i+1,j,n)]+36*u[IND(i+2,j,n)]-16*u[IND(i+3,j,n)]+3*u[IND(i+4,j,n)],0)/(12*dx);
    else if (i+4 > n-1)
        uy = MAX(25*u[IND(i,j,n)]-48*u[IND(i-1,j,n)]+36*u[IND(i-2,j,n)]-16*u[IND(i-3,j,n)]+3*u[IND(i-4,j,n)],0)/(12*dx);
    else
        uy = MAX(MAX(25*u[IND(i,j,n)]-48*u[IND(i+1,j,n)]+36*u[IND(i+2,j,n)]-16*u[IND(i+3,j,n)]+3*u[IND(i+4,j,n)],25*u[IND(i,j,n)]-48*u[IND(i-1,j,n)]+36*u[IND(i-2,j,n)]-16*u[IND(i-3,j,n)]+3*u[IND(i-4,j,n)]),0)/(12*dx);
    rij = sqrt(pow(ux,2)+pow(uy,2))-f[IND(i,j,n)];
    return rij;
}

void Faccurate4thupwindfast(mwSize i, mwSize j, double *u, mwSize n, double *f, double *x, double *y){
    double dx,a,b,c;
    dx = x[1]-x[0];
    if (i-4 < 0)
        a = (48*u[IND(i+1,j,n)]-36*u[IND(i+2,j,n)]+16*u[IND(i+3,j,n)]-3*u[IND(i+4,j,n)])/25;
    else if (i+4 > n-1)
        a = (48*u[IND(i-1,j,n)]-36*u[IND(i-2,j,n)]+16*u[IND(i-3,j,n)]-3*u[IND(i-4,j,n)])/25;
    else
        a = MIN(48*u[IND(i+1,j,n)]-36*u[IND(i+2,j,n)]+16*u[IND(i+3,j,n)]-3*u[IND(i+4,j,n)],48*u[IND(i-1,j,n)]-36*u[IND(i-2,j,n)]+16*u[IND(i-3,j,n)]-3*u[IND(i-4,j,n)])/25;
    if (j-4 < 0)
        b = (48*u[IND(i,j+1,n)]-36*u[IND(i,j+2,n)]+16*u[IND(i,j+3,n)]-3*u[IND(i,j+4,n)])/25;
    else if (j+4 > n-1)
        b = (48*u[IND(i,j-1,n)]-36*u[IND(i,j-2,n)]+16*u[IND(i,j-3,n)]-3*u[IND(i,j-4,n)])/25;
    else
        b = MIN(48*u[IND(i,j+1,n)]-36*u[IND(i,j+2,n)]+16*u[IND(i,j+3,n)]-3*u[IND(i,j+4,n)],48*u[IND(i,j-1,n)]-36*u[IND(i,j-2,n)]+16*u[IND(i,j-3,n)]-3*u[IND(i,j-4,n)])/25;
    c = f[IND(i,j,n)]*12*dx/25;
    u[IND(i,j,n)] = solve(a,b,c);
}

void undivideddifferences(mwSize i, mwSize j, mwSize n, double *u, double *NDD, double *x, mwSize r, double xory){
    mwSize aux,ki,kj;
    if (xory == 1){
        aux = i;
        for (ki=aux-r; ki<aux+r+1; ki++) {
            NDD[IND(ki,0,n)] = u[IND(ki,j,n)];
        }
    }
    else {
        aux = j;
        for (ki=aux-r; ki<aux+r+1; ki++) {
            NDD[IND(ki,0,n)] = u[IND(i,ki,n)];
        }
    }
    for (kj=1;kj<r+1;kj++) {
        for (ki=aux-r;ki<aux+r+1-kj;ki++) {
            NDD[IND(ki,kj,n)] = NDD[IND(ki+1,kj-1,n)]-NDD[IND(ki,kj-1,n)];
        }
    }
}

double findu(mwSize i,mwSize t, mwSize n, mwSize r, double *NDD, double *C, double *x){
    mwSize imin,k;
    double usign,dx;
    dx = x[1]-x[0];
    imin = t;
    for (k=2;k<r+1;k++) {
        if (absolutevalue(NDD[imin+n*k]) >= absolutevalue(NDD[imin-1+n*k]))
            imin = imin-1;   
    }
    usign = 0;
    for (k=1;k<r+1;k++){
        usign = usign + (C[imin-i+r+n*k])*(NDD[imin+n*k]);
    }
    usign = usign/dx;
    return usign;
}
 
double FaccurateENO(mwSize i, mwSize j, double *u, mwSize n, double *f, double *x, double *NDD, double *C, mwSize r){
    double Fa,ux,uplusx,uminusx,uy,uplusy,uminusy,dx;
    dx = x[1]-x[0];
    if (i-r < 0 || i+r > n-1){
        if (r == 2)
            if (i - 2 < 0)
                uy = MAX(3*u[IND(i,j,n)]-4*u[IND(i+1,j,n)]+u[IND(i+2,j,n)],0)/(2*dx);
            else
                uy = MAX(3*u[IND(i,j,n)]-4*u[IND(i-1,j,n)]+u[IND(i-2,j,n)],0)/(2*dx);
        else if (r == 3)
            if (i-3 < 0)
                uy = MAX(11*u[IND(i,j,n)]-18*u[IND(i+1,j,n)]+9*u[IND(i+2,j,n)]+2*u[IND(i+3,j,n)],0)/(6*dx);
            else
                uy = MAX(11*u[IND(i,j,n)]-18*u[IND(i-1,j,n)]+9*u[IND(i-2,j,n)]-2*u[IND(i-3,j,n)],0)/(6*dx);
        else
            if (i-4 < 0)
                uy = (-25*u[IND(i,j,n)]+48*u[IND(i+1,j,n)]-36*u[IND(i+2,j,n)]+16*u[IND(i+3,j,n)]-3*u[IND(i+4,j,n)])/(12*dx);
            else
                uy = (25*u[IND(i,j,n)]-48*u[IND(i-1,j,n)]+36*u[IND(i-2,j,n)]-16*u[IND(i-3,j,n)]+3*u[IND(i-4,j,n)])/(12*dx);
    }
    else{
        undivideddifferences(i,j,n,u,NDD,x,r,1);
        uplusy = findu(i,i,n,r,NDD,C,x);
        uminusy = findu(i,i-1,n,r,NDD,C,x);
        uy = MAX(-uplusy,uminusy);
    }
    if (j-r < 0 || j+r > n-1){
        if (r == 2)
            if (j - 2 < 0)
                ux = MAX(3*u[IND(i,j,n)]-4*u[IND(i,j+1,n)]+u[IND(i,j+2,n)],0)/(2*dx);
            else
                ux = MAX(3*u[IND(i,j,n)]-4*u[IND(i,j-1,n)]+u[IND(i,j-2,n)],0)/(2*dx);
        else if (r == 3)
            if (j-3 < 0)
                ux = MAX(11*u[IND(i,j,n)]-18*u[IND(i,j+1,n)]+9*u[IND(i,j+2,n)]+2*u[IND(i,j+3,n)],0)/(6*dx);
            else
                ux = MAX(11*u[IND(i,j,n)]-18*u[IND(i,j-1,n)]+9*u[IND(i,j-2,n)]-2*u[IND(i,j-3,n)],0)/(6*dx);
        else
            if (j-4 < 0)
                ux = (-25*u[IND(i,j,n)]+48*u[IND(i,j+1,n)]-36*u[IND(i,j+2,n)]+16*u[IND(i,j+3,n)]-3*u[IND(i,j+4,n)])/(12*dx);
            else
                ux = (25*u[IND(i,j,n)]-48*u[IND(i,j-1,n)]+36*u[IND(i,j-2,n)]-16*u[IND(i,j-3,n)]+3*u[IND(i,j-4,n)])/(12*dx);
    }
    else{
        undivideddifferences(i,j,n,u,NDD,x,r,0);
        uplusx = findu(j,j,n,r,NDD,C,x);
        uminusx = findu(j,j-1,n,r,NDD,C,x);
        ux = MAX(-uplusx,uminusx);
    }
    Fa = sqrt(pow(ux,2) + pow(uy,2))-f[IND(i,j,n)];
    return Fa;
}

void FaccurateENOfast(mwSize i, mwSize j, double *u, mwSize n, double *f, double *x, double *NDD, double *C, mwSize r){
    double Fa,ux,uplusx,uminusx,uplusy,uminusy,dx,a,b,c;
    dx = x[1]-x[0];
    if (i-r < 0 || i+r > n-1){
        if (r==2)
            if (i - 2 < 0)
                b = u[IND(i,j,n)]+dx*(-3*u[IND(i,j,n)]+4*u[IND(i+1,j,n)]-u[IND(i+2,j,n)])/(2*dx);
            else
                b = u[IND(i,j,n)]-dx*(3*u[IND(i,j,n)]-4*u[IND(i-1,j,n)]+u[IND(i-2,j,n)])/(2*dx);
        else if (r == 3)
            if (i-3 < 0)
                b = u[IND(i,j,n)]+dx*(-11*u[IND(i,j,n)]+18*u[IND(i+1,j,n)]-9*u[IND(i+2,j,n)]+2*u[IND(i+3,j,n)])/(6*dx);
            else
                b = u[IND(i,j,n)]-dx*(11*u[IND(i,j,n)]-18*u[IND(i-1,j,n)]+9*u[IND(i-2,j,n)]-2*u[IND(i-3,j,n)])/(6*dx);
        else
            if (i-4 < 0)
                b = u[IND(i,j,n)]+dx*(-25*u[IND(i,j,n)]+48*u[IND(i+1,j,n)]-36*u[IND(i+2,j,n)]+16*u[IND(i+3,j,n)]-3*u[IND(i+4,j,n)])/(12*dx);
            else
                b = u[IND(i,j,n)]-dx*(25*u[IND(i,j,n)]-48*u[IND(i-1,j,n)]+36*u[IND(i-2,j,n)]-16*u[IND(i-3,j,n)]+3*u[IND(i-4,j,n)])/(12*dx);
    }
    else{
        undivideddifferences(i,j,n,u,NDD,x,r,1);
        uplusy = findu(i,i,n,r,NDD,C,x);
        uminusy = findu(i,i-1,n,r,NDD,C,x);
        b = MIN(u[IND(i,j,n)]-dx*uminusy,u[IND(i,j,n)]+dx*uplusy);
    }
    if (j-r < 0 || j+r > n-1){
        if (r == 2)
            if (j - 2 < 0)
                a = u[IND(i,j,n)]+dx*(-3*u[IND(i,j,n)]+4*u[IND(i,j+1,n)]-u[IND(i,j+2,n)])/(2*dx);
            else
                a = u[IND(i,j,n)]-dx*(3*u[IND(i,j,n)]-4*u[IND(i,j-1,n)]+u[IND(i,j-2,n)])/(2*dx);
        else if (r == 3)
            if (j-3 < 0)
                a = u[IND(i,j,n)]+dx*(-11*u[IND(i,j,n)]+18*u[IND(i,j+1,n)]-9*u[IND(i,j+2,n)]+2*u[IND(i,j+3,n)])/(6*dx);
            else
                a = u[IND(i,j,n)]-dx*(11*u[IND(i,j,n)]-18*u[IND(i,j-1,n)]+9*u[IND(i,j-2,n)]-2*u[IND(i,j-3,n)])/(6*dx);
        else
            if (j-4 < 0)
                a = u[IND(i,j,n)]+dx*(-25*u[IND(i,j,n)]+48*u[IND(i,j+1,n)]-36*u[IND(i,j+2,n)]+16*u[IND(i,j+3,n)]-3*u[IND(i,j+4,n)])/(12*dx);
            else
                a = u[IND(i,j,n)]-dx*(25*u[IND(i,j,n)]-48*u[IND(i,j-1,n)]+36*u[IND(i,j-2,n)]-16*u[IND(i,j-3,n)]+3*u[IND(i,j-4,n)])/(12*dx);
    }
    else{
        undivideddifferences(i,j,n,u,NDD,x,r,0);
        uplusx = findu(j,j,n,r,NDD,C,x);
        uminusx = findu(j,j-1,n,r,NDD,C,x);
        a = MIN(u[IND(i,j,n)]-dx*uminusx,u[IND(i,j,n)]+dx*uplusx);
    }
    c = f[IND(i,j,n)]*dx;
    u[IND(i,j,n)] = solve(a,b,c);
}


double Sfilter(double x){
    double r;
    if (absolutevalue(x) <= 1)
        r = x;
//     else if (1<= x && x <= 2)
//         r = -x+2;
//     else if (-2 <= x && x <= -1)
//         r = -x-2;
    else
        r = 0;
    return r;
}

double Sfilterd(double x){
    double r;
    if (absolutevalue(x) <= 1)
        r = 1;
//     else if (1<= x && x <= 2)
//         r = -1;
//     else if (-2 <= x && x <= -1)
//         r = -1;
    else
        r = 0;
    return r;
}

void update(mwSize i, mwSize j, double *u, mwSize n, double *f, double *x, double *y, double *p, double *order_accurate, double *fast, double *NDD, double *C, mwSize r){
    double uij,Fm,Fa,uaux,epsilon,dt,dx,ux,uy;
    double uL,uR,uT,uB,uxM,uyM,uxA,uyA;
    mwSize i1,i2,j1,j2;
    if (fast[0]==2){
        dx = x[1]-x[0];
        dt = dx/4;
        epsilon = 2*pow(dx,p[0]);
        
        dx = x[1]-x[0];
        i1 = MAX(i-1,0);
        i2 = MIN(n-1,i+1);
        j1 = MAX(j-1,0);
        j2 = MIN(n-1,j+1);
        uij = u[IND(i,j,n)];
        uT = u[IND(i1,j,n)];
        uB = u[IND(i2,j,n)];
        uL = u[IND(i,j1,n)];
        uR = u[IND(i,j2,n)];
        
        uxM = MAX(MAX(uij-uL,uij-uR),0)/dx;
        uyM = MAX(MAX(uij-uT,uij-uB),0)/dx;
        
        if (order_accurate[0] == 1){
            if (j>0 && j<n-1)
                uxA = absolutevalue(uR-uL)/(2*dx);
            else
                uxA = MAX(MAX(uij-uL,uij-uR),0)/dx;
            if (i>0 && i<n-1)
                uyA = absolutevalue(uB-uT)/(2*dx);
            else
                uyA = MAX(MAX(uij-uT,uij-uB),0)/dx;
        }
        else if (order_accurate[0] == 2){
            if (j < 2)
                uxA = MAX(3*u[IND(i,j,n)]-4*u[IND(i,j+1,n)]+u[IND(i,j+2,n)],0)/(2*dx);
            else if (j > n-3)
                uxA = MAX(3*u[IND(i,j,n)]-4*u[IND(i,j-1,n)]+u[IND(i,j-2,n)],0)/(2*dx);
            else
                uxA = MAX(MAX(3*u[IND(i,j,n)]-4*u[IND(i,j+1,n)]+u[IND(i,j+2,n)],3*u[IND(i,j,n)]-4*u[IND(i,j-1,n)]+u[IND(i,j-2,n)]),0)/(2*dx);
            if (i < 2)
                uyA = MAX(3*u[IND(i,j,n)]-4*u[IND(i+1,j,n)]+u[IND(i+2,j,n)],0)/(2*dx);
            else if (i > n-3)
                uyA = MAX(3*u[IND(i,j,n)]-4*u[IND(i-1,j,n)]+u[IND(i-2,j,n)],0)/(2*dx);
            else
                uyA = MAX(MAX(3*u[IND(i,j,n)]-4*u[IND(i+1,j,n)]+u[IND(i+2,j,n)],3*u[IND(i,j,n)]-4*u[IND(i-1,j,n)]+u[IND(i-2,j,n)]),0)/(2*dx);
        }
        else if (order_accurate[0] == 3){
            if (j-3 < 0)
                uxA = MAX(11*u[IND(i,j,n)]-18*u[IND(i,j+1,n)]+9*u[IND(i,j+2,n)]-2*u[IND(i,j+3,n)],0)/(6*dx);
            else if (j+3 > n-1)
                uxA = MAX(11*u[IND(i,j,n)]-18*u[IND(i,j-1,n)]+9*u[IND(i,j-2,n)]-2*u[IND(i,j-3,n)],0)/(6*dx);
            else
                uxA = MAX(MAX(11*u[IND(i,j,n)]-18*u[IND(i,j+1,n)]+9*u[IND(i,j+2,n)]-2*u[IND(i,j+3,n)],11*u[IND(i,j,n)]-18*u[IND(i,j-1,n)]+9*u[IND(i,j-2,n)]-2*u[IND(i,j-3,n)]),0)/(6*dx);
            if (i-3 < 0)
                uyA = MAX(11*u[IND(i,j,n)]-18*u[IND(i+1,j,n)]+9*u[IND(i+2,j,n)]-2*u[IND(i+3,j,n)],0)/(6*dx);
            else if (i+3 > n-1)
                uyA = MAX(11*u[IND(i,j,n)]-18*u[IND(i-1,j,n)]+9*u[IND(i-2,j,n)]-2*u[IND(i-3,j,n)],0)/(6*dx);
            else
                uyA = MAX(MAX(11*u[IND(i,j,n)]-18*u[IND(i+1,j,n)]+9*u[IND(i+2,j,n)]-2*u[IND(i+3,j,n)],11*u[IND(i,j,n)]-18*u[IND(i-1,j,n)]+9*u[IND(i-2,j,n)]-2*u[IND(i-3,j,n)]),0)/(6*dx);
        }
        else{
            if (j-4 < 0)
                uxA = MAX(25*u[IND(i,j,n)]-48*u[IND(i,j+1,n)]+36*u[IND(i,j+2,n)]-16*u[IND(i,j+3,n)]+3*u[IND(i,j+4,n)],0)/(12*dx);
            else if (j+4 > n-1)
                uxA = MAX(25*u[IND(i,j,n)]-48*u[IND(i,j-1,n)]+36*u[IND(i,j-2,n)]-16*u[IND(i,j-3,n)]+3*u[IND(i,j-4,n)],0)/(12*dx);
            else
                uxA = MAX(MAX(25*u[IND(i,j,n)]-48*u[IND(i,j+1,n)]+36*u[IND(i,j+2,n)]-16*u[IND(i,j+3,n)]+3*u[IND(i,j+4,n)],25*u[IND(i,j,n)]-48*u[IND(i,j-1,n)]+36*u[IND(i,j-2,n)]-16*u[IND(i,j-3,n)]+3*u[IND(i,j-4,n)]),0)/(12*dx);
            if (i-4 < 0)
                uyA = MAX(25*u[IND(i,j,n)]-48*u[IND(i+1,j,n)]+36*u[IND(i+2,j,n)]-16*u[IND(i+3,j,n)]+3*u[IND(i+4,j,n)],0)/(12*dx);
            else if (i+4 > n-1)
                uyA = MAX(25*u[IND(i,j,n)]-48*u[IND(i-1,j,n)]+36*u[IND(i-2,j,n)]-16*u[IND(i-3,j,n)]+3*u[IND(i-4,j,n)],0)/(12*dx);
            else
                uyA = MAX(MAX(25*u[IND(i,j,n)]-48*u[IND(i+1,j,n)]+36*u[IND(i+2,j,n)]-16*u[IND(i+3,j,n)]+3*u[IND(i+4,j,n)],25*u[IND(i,j,n)]-48*u[IND(i-1,j,n)]+36*u[IND(i-2,j,n)]-16*u[IND(i-3,j,n)]+3*u[IND(i-4,j,n)]),0)/(12*dx);
        }
        
        ux = uxM+epsilon*Sfilter((uxA-uxM)/epsilon);
        uy = uyM+epsilon*Sfilter((uyA-uyM)/epsilon);
        uaux = sqrt(pow(ux,2)+pow(uy,2))-f[IND(i,j,n)];
        uij = u[IND(i,j,n)];
        u[IND(i,j,n)] = uij - dt*uaux;
    }
    else {
        dx = x[1]-x[0];
        dt = dx/4;
        epsilon = pow(dx,p[0]);
        Fm = Fmonotone(i,j,u,n,f,x,y);
        if (order_accurate[0] == 1)
            Fa = Faccurate2nd(i,j,u,n,f,x,y);
        else if (order_accurate[0] == 2)
            Fa = Faccurate2ndupwind(i,j,u,n,f,x,y);
        else if (order_accurate[0] == 3)
            Fa = Faccurate3rdupwind(i,j,u,n,f,x,y);
        else if (order_accurate[0] == 4)
            Fa = Faccurate4thupwind(i,j,u,n,f,x,y);
        else {
            Fa = FaccurateENO(i,j,u,n,f,x,NDD,C,r);
        }
        if (fast[0] == 1){
            if (Sfilterd((Fa-Fm)/epsilon) == 1){
                if (order_accurate[0] == 2)
                    Faccurate2ndupwindfast(i,j,u,n,f,x,y);
                else if (order_accurate[0] == 3)
                    Faccurate3rdupwindfast(i,j,u,n,f,x,y);
                else if (order_accurate[0] == 4)
                    Faccurate4thupwindfast(i,j,u,n,f,x,y);
                else
                    FaccurateENOfast(i,j,u,n,f,x,NDD,C,r);
            }
            else
                Fmonotonefast(i,j,u,n,f,x,y);
        }
        else{
            uaux = Fm+epsilon*Sfilter((Fa-Fm)/epsilon);
            uij = u[IND(i,j,n)];
            u[IND(i,j,n)] = uij - dt*uaux;
        }
    }
}

void loopFilter(double *u, double *f, mwSize n, double *dp, mwSize ndp, double *x, double *y, double *p, double *order_accurate, double *fast, double *NDD, double *C,  mwSize r){
    mwSize i,j,k;
    /* loop left, right, up, down */
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            if (dp[IND(i,j,ndp)] == 0){
                update(i,j,u,n,f,x,y,p,order_accurate,fast,NDD,C,r);
            }
        }
    }
    /* loop left, right, down, up */
    for (i=0; i<n; i++) {
        for (j=n-1; j>-1; j--) {
            if (dp[IND(i,j,ndp)] == 0){
                update(i,j,u,n,f,x,y,p,order_accurate,fast,NDD,C,r);
            }
        }
    }
    /* loop right, left, up, down */
    for (i=n-1; i>-1; i--) {
        for (j=0; j<n; j++) {
            if (dp[IND(i,j,ndp)] == 0){
                update(i,j,u,n,f,x,y,p,order_accurate,fast,NDD,C,r);
            }
        }
    }
    
    /* loop right, left, down, up */
    for (i=n-1; i>-1; i--) {
        for (j=n-1; j>-1; j--) {
            if (dp[IND(i,j,ndp)] == 0){
                update(i,j,u,n,f,x,y,p,order_accurate,fast,NDD,C,r);
            }
        }
    }
    
    /* loop up, down, left, right */
    for (j=0; j<n; j++) {
        for (i=0; i<n; i++) {
            if (dp[IND(i,j,ndp)] == 0){
                update(i,j,u,n,f,x,y,p,order_accurate,fast,NDD,C,r);
            }
        }
    }
    /* loop up, down, right, left */
    for (j=0; j<n; j++) {
        for (i=n-1; i>-1; i--) {
            if (dp[IND(i,j,ndp)] == 0){
                update(i,j,u,n,f,x,y,p,order_accurate,fast,NDD,C,r);
            }
        }
    }
    /* loop down, up, left, right */
    for (j=n-1; j>-1; j--) {
        for (i=0; i<n; i++) {
            if (dp[IND(i,j,ndp)] == 0){
                update(i,j,u,n,f,x,y,p,order_accurate,fast,NDD,C,r);
            }
        }
    }
    /* loop down, up, right, left */
    for (j=n-1; j>-1; j--) {
        for (i=n-1; i>-1; i--) {
            if (dp[IND(i,j,ndp)] == 0){
                update(i,j,u,n,f,x,y,p,order_accurate,fast,NDD,C,r);
            }
        }
    }
    
//     for (i=0; i<n; i++) {
//         for (k=0; k<i+1;k++){
//             if (dp[IND(i-k,k,ndp)] == 0){
//                 update(i-k,k,u,n,f,x,y,p,order_accurate,fast,NDD,C,r);
//             }
//         }
//     }
//     for (j=1; j<n; j++){
//         for (k=0; k<n-j; k++) {
//             if (dp[IND(i-n+k,j+k,ndp)] == 0){
//                 update(i-n+k,j+k,u,n,f,x,y,p,order_accurate,fast,NDD,C,r);
//             }
//         }
//     }
}



void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *u, *f, *dp, *x, *y, *unew, *p, *order_accurate, *fast, *NDD, *C, *rv;
  mwSize n,ndp,i,j,count,r;
  
  /* obtain the pointers to the input data */
  u = mxGetPr(prhs[0]);
  f = mxGetPr(prhs[1]);
  dp = mxGetPr(prhs[2]);
  x = mxGetPr(prhs[3]);
  y = mxGetPr(prhs[4]);
  p = mxGetPr(prhs[5]);
  order_accurate = mxGetPr(prhs[6]);
  fast = mxGetPr(prhs[7]);
  NDD = mxGetPr(prhs[8]);
  C = mxGetPr(prhs[9]);
  rv = mxGetPr(prhs[10]);
  
  if (rv[0] == 1)
      r = 1;
  else if (rv[0] == 2)
      r = 2;
  else if (rv[0] == 3)
      r = 3;
  else if (rv[0] == 4)
      r = 4;
  else
      r = 5;
  
  n = mxGetM(prhs[0]);
  ndp = mxGetM(prhs[2]);
  
  /*  call the C subroutine */
  loopFilter(u,f,n,dp,ndp,x,y,p,order_accurate,fast,NDD,C,r);
  
  plhs[0] = mxCreateDoubleMatrix(n,n, mxREAL);
  unew = mxGetPr(plhs[0]);
  /* copy the updated variable u and argmin to the output variable */
  count = 0;
  for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
          unew[count] = u[count];
          NDD[count] = 0;
          count++;
      }
  }
  
}