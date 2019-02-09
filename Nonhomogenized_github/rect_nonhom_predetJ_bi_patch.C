 	#include "mex.h"
	#include "math.h"

	void fun_detJ(double XX[], double YY[], double Z[], double xi[],  double detJ[][1])
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
    double t130, t131, t132, t133, t134, t135, t136, t137, t138, t139;
    double t140, t141, t142, t143, t144, t145, t146, t147, t148, t149;
    double t150, t151, t152, t153, t154, t155, t156, t157, t158, t159;
    double t160, t161, t162, t163, t164, t165, t166, t167, t168, t169;
    double t170, t171, t172, t173, t174, t175, t176, t177, t178, t179;
    double t180, t181, t182, t183, t184, t185, t186, t187, t188, t189;
    double t190, t191, t192, t193, t194, t195, t196, t197, t198, t199;
    double t200, t201, t202, t203, t204, t205, t206, t207, t208, t209;
    double t210, t211, t212, t213, t214, t215, t216, t217, t218, t219;
    double t220, t221, t222, t223, t224, t225, t226, t227, t228, t229;
    double t230, t231, t232, t233, t234, t235, t236, t237, t238, t239;
    double t240, t241, t242, t243, t244, t245, t246, t247, t248, t249;
    double t250, t251, t252, t253, t254, t255, t256, t257, t258, t259;
    double t260, t261, t262, t263, t264, t265, t266, t267, t268, t269;
    double t270, t271, t272, t273, t274, t275, t276, t277, t278, t279;
    double t280, t281, t282, t283, t284, t285, t286, t287, t288, t289;
    double t290, t291, t292, t293, t294, t295, t296, t297, t298, t299;
    double t300, t301, t302, t303, t304, t305, t306, t307, t308, t309;
    double t310, t311, t312, t313, t314, t315, t316, t317, t318, t319;
    double t320, t321, t322, t323, t324, t325, t326, t327, t328, t329;
    double t330, t331, t332, t333, t334, t335, t336, t337, t338, t339;
    double t340, t341, t342, t343, t344, t345, t346, t347, t348, t349;
    double t350, t351, t352, t353, t354, t355, t356, t357, t358, t359;
    double t360, t361, t362, t363, t364, t365, t366, t367, t368, t369;
    double t370, t371, t372, t373, t374, t375, t376, t377, t378, t379;
    double t380, t381, t382, t383, t384, t385, t386, t387, t388, t389;
    double t390, t391, t392, t393, t394, t395, t396, t397, t398, t399;
    double t400, t401, t402, t403, t404, t405, t406, t407, t408, t409;
    double t410, t411, t412, t413, t414, t415, t416, t417, t418, t419;
    double t420, t421, t422, t423, t424, t425, t426, t427, t428, t429;
    double t430, t431, t432, t433, t434, t435, t436, t437, t438, t439;
    double t440, t441, t442, t443, t444, t445, t446, t447, t448, t449;
    double t450, t451, t452, t453, t454, t455, t456, t457, t458, t459;
    double t460, t461, t462, t463, t464, t465, t466, t467, t468, t469;
    double t470, t471, t472, t473, t474, t475, t476, t477, t478, t479;
    double t480, t481, t482, t483, t484, t485, t486, t487, t488, t489;
    double t490, t491, t492, t493, t494, t495, t496, t497, t498, t499;
    double MapleGenVar1, MapleGenVar2, MapleGenVar3, MapleGenVar4;
    double MapleGenVar5, MapleGenVar6, MapleGenVar7;


      t1 = YY[1]*YY[1];
      t4 = YY[3]*YY[3];
      t6 = pow(-Z[0]+Z[1],2.0);
      t8 = XX[0]*XX[0];
      t16 = YY[3]*YY[2];
      t21 = XX[2]*YY[3];
      t23 = XX[3]*YY[2];
      t24 = 2.0*t23;
      t27 = -t4-t6;
      t29 = t16+t6;
      t30 = XX[3]*t29;
      t34 = XX[1]*XX[1];
      t37 = XX[3]*XX[3];
      t39 = YY[0]*YY[0];
      t43 = XX[2]*YY[1];
      t45 = 2.0*t21;
      t46 = 4.0*t23;
      t49 = XX[3]*XX[2];
      t52 = XX[3]*YY[3];
      t55 = -t37-t6;
      t57 = t6*YY[3];
      t61 = YY[2]*YY[2];
      t62 = t61+t6;
      t66 = 2.0*XX[2]*t29;
      t68 = 2.0*t62*XX[3];
      t71 = XX[2]*XX[2];
      t74 = t71*YY[3];
      t75 = XX[2]*t23;
      t76 = YY[2]-YY[3];
      t77 = t76*t6;
      t88 = t8*(-2.0*YY[1]*YY[3]+t1+t4+t6)+XX[0]*(-2.0*(YY[1]-YY[3])*(XX[1]-XX
[3])*YY[0]+2.0*XX[1]*(YY[1]*YY[2]-t16-t6)-2.0*t1*XX[2]+YY[1]*(4.0*t21-t24)+2.0*
XX[2]*t27+2.0*t30)+t39*(-2.0*XX[1]*XX[3]+t34+t37+t6)+YY[0]*(-2.0*t34*YY[2]+XX
[1]*(2.0*t43-t45+t46)+2.0*YY[1]*(-t49-t6)+2.0*XX[2]*t52+2.0*YY[2]*t55+2.0*t57)+
t34*t62+XX[1]*(-2.0*YY[2]*t43+t66-t68)+t1*(t71+t6)+2.0*YY[1]*(-t74+t75+t77)-t71
*t27-2.0*XX[2]*t30-t61*t55-2.0*YY[2]*t57+t6*(t37+t4);
      t89 = xi[1]*xi[1];
      t92 = -2.0*xi[0]*t76;
      t93 = 2.0*YY[2];
      t94 = 2.0*YY[3];
      t100 = 2.0*t16;
      t101 = 2.0*t6;
      t105 = 2.0*xi[0]*t76;
      t108 = XX[2]-XX[3];
      t109 = 2.0*xi[0]*t108;
      t110 = 2.0*XX[2];
      t111 = 2.0*XX[3];
      t120 = 4.0*YY[2];
      t125 = 2.0*t61;
      t126 = 4.0*t6;
      t129 = -2.0*xi[0]*t108;
      t130 = 4.0*XX[2];
      t133 = YY[2]*XX[2];
      t135 = 2.0*xi[0]*(t133-t52);
      t141 = t21-t23;
      t144 = 2.0*xi[0]*t141*t76;
      t153 = 2.0*t49;
      t169 = xi[0]*t141*t108;
      t177 = xi[0]*t76*YY[2];
      t191 = xi[0]*t108*XX[2];
      t199 = t105-t93;
      t201 = t76*t76;
      t202 = xi[0]*xi[0];
      t203 = t202*t201;
      t204 = 2.0*t177;
      t207 = xi[0]*t76;
      t208 = YY[1]+t207-YY[2];
      t209 = xi[0]*t108;
      t220 = t209-XX[2];
      t221 = YY[1]*t220;
      t228 = t108*t108;
      t229 = t202*t228;
      t230 = 2.0*t191;
      t238 = t229-t230+t71+t6;
      MapleGenVar3 = t89*t88;
      MapleGenVar5 = xi[1];
      MapleGenVar7 = t8*(-2.0*t1+YY[1]*(t92+t93+t94)+2.0*xi[0]*t76*YY[3]-t100-
t101)+XX[0]*(YY[0]*(XX[1]*(4.0*YY[1]+t105-t93-t94)+YY[1]*(t109-t110-t111)+2.0*
xi[0]*(-t21-(YY[2]-t94)*XX[3])+t45+t24)+XX[1]*(YY[1]*(t105-t120)+2.0*xi[0]*(-
t61+t4)+t125+t100+t126)+t1*(t129+t130)+YY[1]*(t135+XX[2]*(-t93-4.0*YY[3])+t24)-
t144+t66-t68)+t39*(-2.0*t34+XX[1]*(t129+t110+t111)+2.0*xi[0]*t108*XX[3]-t153-
t101);
      MapleGenVar6 = MapleGenVar7+YY[0]*(t34*(t92+t120)+XX[1]*(YY[1]*(t109-t130
)+t135-2.0*XX[2]*t76-t46)+YY[1]*(2.0*xi[0]*(-t71+t37)+2.0*t71+t153+t126)+2.0*
t169-2.0*t74+2.0*t75+2.0*t77)+2.0*t34*(t177-t61-t6)+XX[1]*(YY[1]*(xi[0]*(XX[2]*
(-t120+t94)+t24)+4.0*t133)+t144-2.0*XX[2]*t29+t68)-2.0*(YY[1]*(-t191+t71+t6)+
t169-t74+t75+t77)*YY[1];
      MapleGenVar4 = MapleGenVar5*MapleGenVar6;
      MapleGenVar2 = MapleGenVar3+MapleGenVar4;
      MapleGenVar1 = MapleGenVar2+t8*(YY[1]*t199+t1+t203-t204+t6+t61)+XX[0]*(
-2.0*YY[0]*(XX[1]+t209-XX[2])*t208+XX[1]*(YY[1]*(t92+t93)-2.0*t203+4.0*t177-
t125-t101)+2.0*t208*t221);
      t250 = MapleGenVar1+t39*(t34+XX[1]*(t109-t110)+t229-t230+t71+t6)+YY[0]*(
t34*t199+2.0*XX[1]*(-YY[1]+t207-YY[2])*t220-2.0*t238*YY[1])+t34*(t203-t204+t61+
t6)-2.0*XX[1]*(t207-YY[2])*t221+t238*t1;

   }


	/* The gateway routine */
    void mexFunction(int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray*prhs[])

    {
        double *XX, *YY, *Z, *xi, *Pre_detJ_bi;

		/* Check for proper number of arguments */

		if (nrhs != 4) {
			mexErrMsgTxt("Four input arguments required.");
		}
		else if (nlhs != 1) {
			mexErrMsgTxt("One output argument required.");
		}

		XX = mxGetPr(prhs[0]);
        YY = mxGetPr(prhs[1]);
		Z = mxGetPr(prhs[2]);
		xi = mxGetPr(prhs[3]);


		/* Create a matrix for the return argument */
		plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);

		/* Assign pointers to each input and output. */
		Pre_detJ_bi = mxGetPr(plhs[0]);

		/* Call the C subroutine */
        fun_detJ(XX,YY,Z,xi,Pre_detJ_bi);
	}
