 	#include "mex.h"
	#include "math.h"

	void prisma_M(double X[], double Y[], double Z[], double lambda[],  double M[][6])
	{
	
    double t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;
    double t10, t11, t12, t13, t14, t15, t16, t17, t18, t19;
    double t20, t21, t22, t23, t24, t25, t26, t27, t28, t29;
    double t30, t31, t32, t33, t34, t35, t36, t37, t38, t39;
    double t40, t41, t42, t43, t44, t45, t46, t47, t48, t49;
    double t50, t51, t52, t53, t54, t55, t56, t57, t58, t59;
    double t60, t61, t62, t63, t64, t65, t66, t67, t68, t69;
    double t70, t71, t72, t73, t74, t75, t76, t77, t78, t79;
    double t80, t81, t82, t83, t84, t85, t86, t87, t88, t89;
    double t90, t91, t92, t93, t94, t95, t96, t97, t98, t99;
    
    
    
	  
      t1 = Z[0]-Z[1];
      t2 = lambda[0]*lambda[0];
      t4 = lambda[1]*lambda[0];
      t5 = 3.0*t4;
      t6 = lambda[1]*lambda[1];
      t15 = X[0]*Y[1]-X[0]*Y[2]-X[1]*Y[0]+X[1]*Y[2]+X[2]*Y[0]-X[2]*Y[1];
      t16 = t15*(6.0*t2+t5+t6)*t1;
      t17 = t16/360.0;
      t18 = t16/720.0;
      t24 = t15*(3.0*t2+4.0*t4+3.0*t6)*t1;
      t25 = t24/720.0;
      t26 = t24/1440.0;
      t30 = t15*(t2+t5+6.0*t6)*t1;
      t31 = t30/360.0;
      t32 = t30/720.0;
      M[0][0] = -t17;
      M[0][1] = -t18;
      M[0][2] = -t18;
      M[0][3] = -t25;
      M[0][4] = -t26;
      M[0][5] = -t26;
      M[1][0] = -t18;
      M[1][1] = -t17;
      M[1][2] = -t18;
      M[1][3] = -t26;
      M[1][4] = -t25;
      M[1][5] = -t26;
      M[2][0] = -t18;
      M[2][1] = -t18;
      M[2][2] = -t17;
      M[2][3] = -t26;
      M[2][4] = -t26;
      M[2][5] = -t25;
      M[3][0] = -t25;
      M[3][1] = -t26;
      M[3][2] = -t26;
      M[3][3] = -t31;
      M[3][4] = -t32;
      M[3][5] = -t32;
      M[4][0] = -t26;
      M[4][1] = -t25;
      M[4][2] = -t26;
      M[4][3] = -t32;
      M[4][4] = -t31;
      M[4][5] = -t32;
      M[5][0] = -t26;
      M[5][1] = -t26;
      M[5][2] = -t25;
      M[5][3] = -t32;
      M[5][4] = -t32;
      M[5][5] = -t31;
    }


	/* The gateway routine */
    void mexFunction(int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray*prhs[])

    {
        double *X, *Y, *Z, *lambda, *M;
   
		/* Check for proper number of arguments */

		if (nrhs != 4) {
			mexErrMsgTxt("Four input arguments required.");
		}
		else if (nlhs != 1) {
			mexErrMsgTxt("One output argument required.");
		}

		X       = mxGetPr(prhs[0]);
        Y       = mxGetPr(prhs[1]);
		Z       = mxGetPr(prhs[2]);		
		lambda  = mxGetPr(prhs[3]);


		
		/* Create a matrix for the return argument */
		plhs[0] = mxCreateDoubleMatrix(6,6,mxREAL);

		/* Assign pointers to each input and output. */
		M = mxGetPr(plhs[0]);

		/* Call the C subroutine */
        prisma_M(X,Y,Z,lambda,M);
	}    
    
    
    


