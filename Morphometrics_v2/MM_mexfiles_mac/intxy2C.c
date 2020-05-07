#include <math.h>
//#include <iostream.h>
#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    double *ax,*ay,*bx,*by,*jstartp,*jdirp;
    double *ai,*bi,*xi,*yi;
    mwSize m,n,k;
    int i,j,cl,d,dtmp,jstart,jdir,p;
    double ux,uy,vx,vy,wx,wy,vvx,vvy,zx,zy,cpuv,cpuvv,cpvz,cpwz,*al,*bl,*xl,*yl,intx,inty,ia,ib,magu,magv,magz,r;
    bool f,o,multiout;

    /* Check for proper number of arguments */
    if (nrhs != 4) { 
	mexErrMsgTxt("4 input arguments required."); 
    } else if (nlhs != 1) {
	mexErrMsgTxt("1 output argument required."); 
    } 
    /* Check the dimensions of Y.  Y can be 4 X 1 or 1 X 4. */ 
    
    m = mxGetN(prhs[0]); 
    n = mxGetM(prhs[0]);
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || !mxIsDouble(prhs[0])
    || mxIsComplex(prhs[1]) || ((n != 1) && (m != 1))) { 
	mexErrMsgTxt("Inputs 1-4 must be Nx1 real vectors."); 
    }
    if(n>m)m=n;
    m = m/2;
    k=mxGetN(prhs[2]);
    if(mxGetM(prhs[2])>k)k=mxGetM(prhs[2]);
    
    /* Assign pointers to the various parameters */ 
    ax = mxGetPr(prhs[0]);
    ay = mxGetPr(prhs[1]);
    bx = mxGetPr(prhs[2]);
    by = mxGetPr(prhs[3]);

    /* Create a matrix for the return argument */ 
    plhs[0] = mxCreateDoubleMatrix(1, m, mxREAL);
    ai = mxGetPr(plhs[0]);
    
    /* Do the actual computations */
    for(i=0;i<m;i++)
    {
        ai[i]=0;
        ux = ax[2*i+1]-ax[2*i];
        uy = ay[2*i+1]-ay[2*i];
        for(j=0;j<k;j++)
        {
            if((ax[2*i]==bx[j])&&(ay[2*i]==by[j])) 
            {
                ai[i]=1;
                continue;
            }
            if(j>0)
            {
                vx = bx[j]-ax[2*i];
                vy = by[j]-ay[2*i];
                wx = bx[j]-ax[2*i+1];
                wy = by[j]-ay[2*i+1];
                vvx = bx[j-1]-ax[2*i];
                vvy = by[j-1]-ay[2*i];
                zx = bx[j]-bx[j-1];
                zy = by[j]-by[j-1];
                //printf("cl=%d, i=%d, j=%d, k=%d\n",cl,i,j,k);
                cpuv = ux*vy-uy*vx;
                cpuvv = ux*vvy-uy*vvx;
                cpvz = vx*zy-vy*zx;
                cpwz = wx*zy-wy*zx;
                if((cpuv*cpuvv<0)&&(cpvz*cpwz<0))
                {
                    ai[i]=1;
                }
            }
        }
    }
    return;
}