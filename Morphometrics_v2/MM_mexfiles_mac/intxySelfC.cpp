#include <math.h>
#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    double *ax,*ay; 
    double *ai,*bi;
    mwSize m,n;
    int i,j,*al,*bl,cl,d,dtmp,maxnint;
    double ux,uy,vx,vy,wx,wy,vvx,vvy,zx,zy,cpuv,cpuvv,cpvz,cpwz;

    /* Check for proper number of arguments */
    
    if (nrhs != 2) { 
	mexErrMsgTxt("Two input arguments required."); 
    } else if (nlhs != 2) {
	mexErrMsgTxt("Two output arguments required."); 
    } 
    
    /* Check the dimensions of Y.  Y can be 4 X 1 or 1 X 4. */ 
    
    m = mxGetM(prhs[0]); 
    n = mxGetN(prhs[0]);
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || !mxIsDouble(prhs[0])
    || mxIsComplex(prhs[1]) || (n != 1)) { 
	mexErrMsgTxt("Inputs must be Nx1 real vectors."); 
    } 
    
    /* Assign pointers to the various parameters */ 
   
    ax = mxGetPr(prhs[0]); 
    ay = mxGetPr(prhs[1]);
    /* Do the actual computations in a subroutine */
    maxnint = 2*m;
    al = new int [maxnint];
    bl = new int [maxnint];
    cl = 0;
    for(i=0;i<m;i++)
    {
        ux = ax[(i+1)%m]-ax[i];
        uy = ay[(i+1)%m]-ay[i];
        for(j=0;j<m;j++)
        {
            if(i==j) continue;
            if((ax[i]==ax[j])&&(ay[i]==ay[j])) 
            {
                al[cl]=i;
                bl[cl]=j;
                cl++;
                continue;
            }
            vx = ax[j]-ax[i];
            vy = ay[j]-ay[i];
            wx = ax[j]-ax[(i+1)%m];
            wy = ay[j]-ay[(i+1)%m];
            vvx = ax[(j+m-1)%m]-ax[i];
            vvy = ay[(j+m-1)%m]-ay[i];
            zx = ax[j]-ax[(j+m-1)%m];
            zy = ay[j]-ay[(j+m-1)%m];
            cpuv = ux*vy-uy*vx;
            cpuvv = ux*vvy-uy*vvx;
            cpvz = vx*zy-vy*zx;
            cpwz = wx*zy-wy*zx;
            if((cpuv*cpuvv<0)&&(cpvz*cpwz<0))
            {
                if(cl==maxnint){i=m;cl=0;break;}
                al[cl]=i;
                bl[cl]=j;
                cl++;
            }
        }
    }
    
    /* Create a matrix for the return argument */ 
    plhs[0] = mxCreateDoubleMatrix(cl, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(cl, 1, mxREAL);
    ai = mxGetPr(plhs[0]);
    bi = mxGetPr(plhs[1]);
    /*for(i=0;i<cl;i++)
    {
        ai[i] = al[i]+1;
        bi[i] = bl[i]+1;
    }*/
    if(cl==0){delete[] al; delete[] bl; return;}
    d = 0;
    for(i=0;i<cl;i++)
    {
        dtmp = (al[i]-al[(i+cl-1)%cl]+m)%m;
        //printf("i=%d, dtmp=%d, al[i]=%d, bl[i]=%d\n",i,dtmp,al[i],bl[i]);
        if(dtmp>d)
        {
            d = dtmp;
            j = i;
        }
    }
    //printf("cl=%d, d=%d, j=%d\n",cl,d,j);
    for(i=0;i<cl;i++)
    {
        ai[i] = al[(i-j+cl)%cl]+1;
        bi[i] = (bl[(i-j+cl)%cl]+m-1)%m+1;//bl[(i-j+cl)%cl]+1;//
        //printf("al[i]=%d, ai[i]=%f, bl[i]=%d, bi[i]=%f, i=%d, m=%d, espr=%d, espr=%d\n",al[i],ai[i],bl[i],bi[i],i,m,bl[(i-j+cl)%cl],(bl[(i-j+cl)%cl]+m-1)%m+1);
    }
    delete[] al;
    delete[] bl;
    return;
}