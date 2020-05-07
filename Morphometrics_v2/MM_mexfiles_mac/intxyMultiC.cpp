#include <math.h>
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
    if (nrhs != 6) { 
	mexErrMsgTxt("6 input arguments required."); 
    } else if (nlhs != 4) {
	mexErrMsgTxt("4 output arguments required."); 
    } 
    /* Check the dimensions of Y.  Y can be 4 X 1 or 1 X 4. */ 
    
    m = mxGetN(prhs[0]); 
    n = mxGetM(prhs[0]);
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || !mxIsDouble(prhs[0])
    || mxIsComplex(prhs[1]) || ((n != 1) && (m != 1))) { 
	mexErrMsgTxt("Inputs 1-4 must be Nx1 real vectors."); 
    }
    if(n>m)m=n;
    k=mxGetN(prhs[2]);
    if(mxGetM(prhs[2])>k)k=mxGetM(prhs[2]);
    
    /* Assign pointers to the various parameters */ 
    ax = mxGetPr(prhs[0]);
    ay = mxGetPr(prhs[1]);
    bx = mxGetPr(prhs[2]);
    by = mxGetPr(prhs[3]);
    jstartp = mxGetPr(prhs[4]);
    jdirp = mxGetPr(prhs[5]);
    jstart = *jstartp-1;
    jdir = *jdirp;
    multiout = false;
    if(jstart<0){multiout=true;jstart=0;}
    /* Do the actual computations in a subroutine */
    al = new double [m*k];
    bl = new double [m*k];
    xl = new double [m*k];
    yl = new double [m*k];
    cl = 0;
    o = true;
    for(i=0;i<m-1;i++)
    {
        ux = ax[(i+1)%m]-ax[i];
        uy = ay[(i+1)%m]-ay[i];
        magu=sqrt(ux*ux+uy*uy);
        if(magu==0)continue;
        ux=ux/magu;
        uy=uy/magu;
        f=true;
        p=jstart;
        //printf("p=%d, jstart=%d, jdir=%d, k=%d, test=%d\n",p,jstart,jdir,k,-1%2);
        while((p!=jstart)|f)
        {
            j=p;
            f=false;
            p=(p+jdir+k)%k;
            if((ax[i]==bx[j])&&(ay[i]==by[j])) 
            {
                al[cl]=i;
                bl[cl]=j;
                xl[cl]=ax[i];
                yl[cl]=ay[i];
                cl++;
                continue;
            }
            if(j>0)
            {
                vx = bx[j]-ax[i];
                vy = by[j]-ay[i];
                wx = bx[j]-ax[i+1];
                wy = by[j]-ay[i+1];
                vvx = bx[j-1]-ax[i];
                vvy = by[j-1]-ay[i];
                zx = bx[j]-bx[j-1];
                zy = by[j]-by[j-1];
                //printf("cl=%d, i=%d, j=%d, k=%d\n",cl,i,j,k);
                //printf("i=%d, j=%d, vx=%f, vy=%f, wx=%f, wy=%f, vvx=%f, vvy=%f\n",i,j,vx,vy,wx,wy,vvx,vvy);
                cpuv = ux*vy-uy*vx;
                cpuvv = ux*vvy-uy*vvx;
                cpvz = vx*zy-vy*zx;
                cpwz = wx*zy-wy*zx;
                if((cpuv*cpuvv<=0)&&(cpvz*cpwz<0))
                {
                    magv=sqrt(vx*vx+vy*vy);
                    if(magv==0)continue;
                    vx=vx/magv;
                    vy=vy/magv;
                    magz=sqrt(zx*zx+zy*zy);
                    if(magz==0)continue;
                    zx=zx/magz;
                    zy=zy/magz;
                    r=magv*sin(acos(vx*zx+vy*zy))/sin(acos(-ux*zx-uy*zy));
                    if((cpuv*cpuvv==0)&&(r==1))continue;
                    ia=i+r/magu;
                    intx=ax[i]+r/magu*(ax[i+1]-ax[i]);
                    inty=ay[i]+r/magu*(ay[i+1]-ay[i]);
                    r=sqrt((bx[j]-intx)*(bx[j]-intx)+(by[j]-inty)*(by[j]-inty));
                    if((cpuv*cpuvv==0)&&(r==1))continue;
                    ib=j-r/magz;
                    al[cl]=ia;
                    bl[cl]=ib;
                    xl[cl]=intx;
                    yl[cl]=inty;
                    cl++;
                    //printf("i=%d, j=%d, magv=%f, magz=%f, magu=%f, r=%f, zx=%f, zy=%f\n",i,j,magv,magz,magu,r,zx,zy);
                    if(!multiout)
                    {
                        o=false;
                        break;
                    }
                }
            }
        }
        if(!o)break;
    }
    
    /* Create a matrix for the return argument */ 
    plhs[0] = mxCreateDoubleMatrix(1, cl, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, cl, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, cl, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1, cl, mxREAL);
    ai = mxGetPr(plhs[2]);
    bi = mxGetPr(plhs[3]);
    xi = mxGetPr(plhs[0]);
    yi = mxGetPr(plhs[1]);
    //printf("cl=%d, i=%d, j=%d\n",cl,i,j);
    if(cl==0)return;
    // Now asign the outputs
    for(i=0;i<cl;i++)
    {
        ai[i] = al[i]+1;
        bi[i] = bl[i]+1;
        xi[i] = xl[i];
        yi[i] = yl[i];
    }
    delete [] al;
    delete [] bl;
    delete [] xl;
    delete [] yl;
    return;
}