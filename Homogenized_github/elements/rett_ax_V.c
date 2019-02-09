 	#include "mex.h"
	#include "math.h"

	void rett_V(double S[], double Z[], double V[])
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
	  
	  
      t3 = S[0]-S[1];
      t5 = Z[0]-Z[1];
      t7 = (2.0*S[0]+S[1])*t3*t5/12.0;
      t12 = (S[0]+2.0*S[1])*t3*t5/12.0;
      V[0] = t7;
      V[1] = t12;
      V[2] = t7;
      V[3] = t12;
    }


	/* The gateway routine */
    void mexFunction(int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray*prhs[])

    {
        double *S, *Z, *V;
   
		/* Check for proper number of arguments */

		if (nrhs != 2) {
			mexErrMsgTxt("Two input arguments required.");
		}
		else if (nlhs != 1) {
			mexErrMsgTxt("One output argument required.");
		}

		S = mxGetPr(prhs[0]);
        Z = mxGetPr(prhs[1]);

		
		/* Create a matrix for the return argument */
		plhs[0] = mxCreateDoubleMatrix(4,1,mxREAL);

		/* Assign pointers to each input and output. */
		V = mxGetPr(plhs[0]);

		/* Call the C subroutine */
	    rett_V(S, Z, V);
	}
