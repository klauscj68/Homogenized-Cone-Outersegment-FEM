 	#include "mex.h"
	#include "math.h"

	void rectangule_M_gen( double X[], double Y[], double Z[], double lambda[],  double M[][4])
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
    double t100, t101, t102, t103, t104, t105, t106, t107, t108, t109;
    double t110, t111, t112, t113, t114, t115, t116, t117, t118, t119;
    double t120, t121, t122, t123, t124, t125, t126, t127, t128, t129;
	  
      t3 = X[0]*X[0];
      t4 = Y[1]*Y[1];
      t5 = t4*t3;
      t6 = lambda[0]*lambda[0];
      t8 = lambda[0]*lambda[1];
      t11 = lambda[1]*lambda[1];
      t13 = X[0]*X[1];
      t14 = Y[0]*Y[1];
      t26 = X[1]*X[1];
      t27 = Y[0]*Y[0];
      t28 = t27*t26;
      t33 = Z[0]*Z[0];
      t38 = Z[1]*Z[1];
      t42 = 4.0*Y[1]*lambda[0]*lambda[1]*Y[0]*t13-2.0*t11*t14*t13-2.0*t6*t14*
t13-2.0*Z[1]*Z[0]*t3+t11*t28+t11*t5-2.0*t33*t13+t6*t28-2.0*t8*t28+t33*t3+t38*t3
+t6*t5-2.0*t5*t8;
      t43 = Z[0]*Z[1];
      t69 = -2.0*Z[1]*Z[0]*t26-2.0*Z[1]*Z[0]*t27-2.0*Z[1]*Z[0]*t4-2.0*t38*t13+
4.0*t43*t13-2.0*t33*t14-2.0*t38*t14+4.0*t43*t14+t33*t26+t38*t26+t33*t27+t38*t27
+t33*t4+t38*t4;
      t71 = sqrt(t42+t69);
      t74 = (3.0*lambda[0]+lambda[1])*t71;
      t75 = t74/36.0;
      t76 = t74/72.0;
      t78 = (lambda[0]+lambda[1])*t71;
      t79 = t78/36.0;
      t80 = t78/72.0;
      t83 = (lambda[0]+3.0*lambda[1])*t71;
      t84 = t83/36.0;
      t85 = t83/72.0;
      M[0][0] = t75;
      M[0][1] = t76;
      M[0][2] = t79;
      M[0][3] = t80;
      M[1][0] = t76;
      M[1][1] = t75;
      M[1][2] = t80;
      M[1][3] = t79;
      M[2][0] = t79;
      M[2][1] = t80;
      M[2][2] = t84;
      M[2][3] = t85;
      M[3][0] = t80;
      M[3][1] = t79;
      M[3][2] = t85;
      M[3][3] = t84;



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
		plhs[0] = mxCreateDoubleMatrix(4,4,mxREAL);

		/* Assign pointers to each input and output. */
		M = mxGetPr(plhs[0]);

		/* Call the C subroutine */
        rectangule_M_gen(X,Y,Z,lambda,M);
	}
