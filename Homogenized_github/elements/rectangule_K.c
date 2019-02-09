 	#include "mex.h"
	#include "math.h"

	void rectangule_K( double L[], double H[], double lambda[],  double K[][4])
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
	  
      t1 = lambda[1]*lambda[1];
      t2 = log(lambda[0]);
      t3 = t2*t1;
      t4 = lambda[0]*lambda[0];
      t5 = H[0]*H[0];
      t6 = t5*t4;
      t8 = 24.0*t6*t3;
      t9 = log(lambda[1]);
      t10 = t9*t1;
      t12 = 24.0*t6*t10;
      t13 = L[0]*L[0];
      t14 = t13*t4;
      t15 = t14*t3;
      t16 = 2.0*t15;
      t17 = t1*lambda[1];
      t19 = lambda[0]*t13;
      t20 = t19*t2*t17;
      t21 = 4.0*t20;
      t22 = t1*t1;
      t25 = 2.0*t13*t2*t22;
      t26 = t14*t10;
      t27 = 2.0*t26;
      t29 = t19*t9*t17;
      t30 = 4.0*t29;
      t33 = 2.0*t13*t9*t22;
      t34 = t4*t4;
      t35 = t34*t5;
      t36 = 12.0*t35;
      t37 = t4*lambda[0];
      t40 = 48.0*lambda[1]*t37*t5;
      t41 = t1*t6;
      t42 = 36.0*t41;
      t44 = 3.0*t34*t13;
      t47 = 6.0*lambda[1]*t37*t13;
      t49 = 6.0*t17*t19;
      t51 = 3.0*t22*t13;
      t52 = t8-t12+t16-t21+t25-t27+t30-t33+t36-t40+t42+t44-t47+t49-t51;
      t53 = 1/L[0];
      t58 = lambda[0]-lambda[1];
      t59 = t58*t58;
      t62 = 1/t59/t58/lambda[0]/H[0];
      t64 = t62*t53*t52/24.0;
      t65 = t8-t12+t16-t21+t25-t27+t30-t33+t36-t40+t42-t44+t47-t49+t51;
      t68 = t62*t53*t65/24.0;
      t69 = t2*t37;
      t70 = t5*lambda[1];
      t72 = 24.0*t70*t69;
      t73 = t9*t37;
      t75 = 24.0*t70*t73;
      t76 = lambda[1]*t13;
      t77 = t76*t69;
      t78 = 2.0*t77;
      t79 = 4.0*t15;
      t80 = 2.0*t20;
      t81 = t76*t73;
      t82 = 2.0*t81;
      t83 = 4.0*t26;
      t84 = 2.0*t29;
      t85 = 12.0*t41;
      t86 = t72-t75+t78-t79+t80-t82+t83-t84-t36+t85+t44-t47+t49-t51;
      t89 = t62*t53*t86/24.0;
      t90 = t72-t75+t78-t79+t80-t82+t83-t84-t36+t85-t44+t47-t49+t51;
      t93 = t62*t53*t90/24.0;
      t94 = t2*t34;
      t96 = 24.0*t5*t94;
      t97 = t9*t34;
      t99 = 24.0*t5*t97;
      t101 = 2.0*t13*t94;
      t102 = 4.0*t77;
      t104 = 2.0*t13*t97;
      t105 = 4.0*t81;
      t106 = 36.0*t35;
      t107 = t96-t99+t101-t102+t16-t104+t105-t27-t106+t40-t85+t44-t47+t49-t51;
      t110 = t62*t53*t107/24.0;
      t111 = t96-t99+t101-t102+t16-t104+t105-t27-t106+t40-t85-t44+t47-t49+t51;
      t114 = t62*t53*t111/24.0;
      K[0][0] = t64;
      K[0][1] = -t68;
      K[0][2] = -t89;
      K[0][3] = t93;
      K[1][0] = -t68;
      K[1][1] = t64;
      K[1][2] = t93;
      K[1][3] = -t89;
      K[2][0] = -t89;
      K[2][1] = t93;
      K[2][2] = t110;
      K[2][3] = -t114;
      K[3][0] = t93;
      K[3][1] = -t89;
      K[3][2] = -t114;
      K[3][3] = t110;

    }


	/* The gateway routine */
    void mexFunction(int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray*prhs[])

    {
        double *L, *H, *lambda, *K;
   
		/* Check for proper number of arguments */

		if (nrhs != 3) {
			mexErrMsgTxt("Three input arguments required.");
		}
		else if (nlhs != 1) {
			mexErrMsgTxt("One output argument required.");
		}

		L       = mxGetPr(prhs[0]);
        H       = mxGetPr(prhs[1]);
		lambda  = mxGetPr(prhs[2]);


		
		/* Create a matrix for the return argument */
		plhs[0] = mxCreateDoubleMatrix(4,4,mxREAL);

		/* Assign pointers to each input and output. */
		K = mxGetPr(plhs[0]);

		/* Call the C subroutine */
        rectangule_K(L,H,lambda,K);
	}
