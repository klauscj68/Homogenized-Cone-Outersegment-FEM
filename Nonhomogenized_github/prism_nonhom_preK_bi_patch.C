 	#include "mex.h"
	#include "math.h"

	void fun_K_nabla_x_u(double XX[], double YY[], double Z[], double xi[],  double Pre_K_output[][18])
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


      t1 = -YY[1]+YY[2]+YY[4]-YY[5];
      t3 = xi[2]*t1+YY[1]-YY[2];
      t4 = xi[2]-1.0;
      t6 = XX[3]-XX[5]-XX[0]+XX[2];
      t7 = YY[1]*t6;
      t8 = XX[0]-XX[1]-XX[3]+XX[4];
      t9 = YY[2]*t8;
      t10 = YY[4]-YY[5];
      t11 = XX[0]*t10;
      t16 = -XX[4]+XX[5];
      t17 = YY[0]*t16;
      t18 = -XX[3]+XX[5];
      t19 = YY[4]*t18;
      t20 = XX[3]-XX[4];
      t21 = YY[5]*t20;
      t22 = -t16*YY[3];
      t24 = xi[2]*xi[2];
      t26 = 2.0*XX[0];
      t27 = 2.0*XX[2];
      t29 = YY[1]*(-XX[3]+XX[5]+t26-t27);
      t30 = 2.0*XX[1];
      t32 = YY[2]*(XX[3]-XX[4]-t26+t30);
      t34 = 2.0*YY[0];
      t39 = -YY[0]*t16;
      t42 = -XX[0]+XX[2];
      t43 = YY[1]*t42;
      t45 = YY[2]*(XX[0]-XX[1]);
      t46 = XX[1]-XX[2];
      t47 = YY[0]*t46;
      t49 = 1/(t24*(t7+t9+t11+(-YY[3]+YY[5]+YY[0])*XX[1]+(YY[3]-YY[4]-YY[0])*XX
[2]+t17+t19+t21+t22)+xi[2]*(t29+t32-XX[0]*t10+XX[1]*(YY[3]-YY[5]-t34)+XX[2]*(-
YY[3]+YY[4]+t34)+t39)+t43+t45+t47);
      t51 = XX[1]-XX[2]-XX[4]+XX[5];
      t53 = xi[2]*t51-XX[1]+XX[2];
      t55 = YY[0]-YY[2]-YY[3]+YY[5];
      t56 = XX[1]*t55;
      t57 = -YY[0]+YY[1]+YY[3]-YY[4];
      t58 = XX[2]*t57;
      t59 = XX[0]*t1;
      t60 = -YY[1]*t18;
      t61 = -YY[2]*t20;
      t62 = YY[3]-YY[5];
      t63 = XX[4]*t62;
      t64 = -YY[3]+YY[4];
      t65 = XX[5]*t64;
      t66 = XX[3]*t10;
      t69 = 2.0*YY[2];
      t70 = -t34+t69+YY[3]-YY[5];
      t72 = 2.0*YY[1];
      t73 = t34-t72-YY[3]+YY[4];
      t75 = t72-t69-YY[4]+YY[5];
      t77 = YY[1]*t18;
      t78 = YY[2]*t20;
      t80 = xi[2]*(XX[1]*t70+XX[2]*t73+XX[0]*t75+t39+t77+t78);
      t81 = YY[0]-YY[2];
      t82 = XX[1]*t81;
      t83 = YY[1]-YY[0];
      t84 = t83*XX[2];
      t85 = YY[1]-YY[2];
      t86 = XX[0]*t85;
      t88 = 1/(t24*(t56+t58+t59+t17+t60+t61+t63+t65-t66)+t80+t82+t84-t86);
      t90 = YY[2]-YY[5];
      t92 = -YY[1]+YY[4];
      t94 = YY[2]*XX[4];
      t95 = XX[4]*YY[5];
      t96 = XX[5]*YY[1];
      t97 = YY[4]*XX[5];
      t100 = 1.0-xi[0]-xi[1];
      t101 = YY[3]*t100;
      t102 = YY[4]*xi[0];
      t103 = YY[5]*xi[1];
      t106 = -YY[3]*t100;
      t109 = -XX[3]*t100;
      t110 = xi[0]*XX[4];
      t111 = xi[1]+1.0;
      t112 = t111*XX[5];
      t115 = XX[3]*t100;
      t116 = xi[1]*XX[5];
      t124 = -t95+t97;
      t132 = t109-t110-t116;
      t141 = t16*YY[3];
      t159 = -t83*XX[2];
      t164 = 1/(Z[0]-Z[1]);
      t167 = -xi[2]*t55+YY[0]-YY[2];
      t169 = YY[0]*t51;
      t170 = YY[4]-YY[5]-YY[1];
      t172 = -XX[1]*t62;
      t173 = YY[3]-YY[4]+YY[1];
      t178 = XX[4]-XX[5]-t30+t27;
      t179 = YY[0]*t178;
      t182 = XX[1]*t62;
      t183 = -YY[3]+YY[4]-t72;
      t189 = 1/(t24*(XX[0]*t170+XX[2]*t173+YY[4]*t18+t169+t172+t21+t22+t60+t9)+
xi[2]*(t179+t32+XX[0]*(-YY[4]+YY[5]+t72)+t182+XX[2]*t183-t60)+t47+t45+YY[1]*t42
);
      t192 = -xi[2]*t6-XX[0]+XX[2];
      t194 = -XX[3]*t10;
      t197 = -XX[0]*t85;
      t199 = 1/(t24*(t59+t58+t56+t17+t60+t61+t194+t65+t63)+t80+t197+t84+t82);
      t202 = XX[2]-XX[5];
      t204 = YY[2]*XX[3];
      t205 = YY[3]*XX[2];
      t206 = YY[3]*XX[5];
      t207 = YY[5]*XX[3];
      t215 = -xi[0]-xi[1]+2.0;
      t239 = -XX[0]*t1;
      t240 = -YY[0]*t51;
      t241 = -XX[2]*t173;
      t243 = (XX[1]+XX[3]-XX[4])*YY[2];
      t248 = XX[1]-XX[4];
      t252 = -XX[0]*t75;
      t253 = -YY[0]*t178;
      t254 = -XX[2]*t183;
      t256 = YY[2]*(-t30-XX[3]+XX[4]);
      t257 = YY[1]*XX[3];
      t258 = YY[3]*XX[1];
      t262 = -YY[0]*t46;
      t263 = YY[1]*XX[2];
      t264 = YY[2]*XX[1];
      t270 = xi[2]*t57+YY[0]-YY[1];
      t272 = YY[4]-YY[5]+YY[2];
      t274 = -YY[2]-YY[3]+YY[5];
      t276 = -XX[2]*t64;
      t281 = YY[3]-YY[5]+t69;
      t283 = XX[2]*t64;
      t287 = 1/(t24*(XX[0]*t272+XX[1]*t274+t169+t19+t21+t22+t276+t61+t7)+xi[2]*
(t179+t29+XX[0]*(-YY[4]+YY[5]-t69)+XX[1]*t281+t283+t78)+t43+t45+t47);
      t290 = xi[2]*t8-XX[0]+XX[1];
      t292 = -XX[5]*t64;
      t296 = 1/(t24*(t59+t56+t58+t17+t60+t61+t194+t63-t292)+t80+t197+t82-t159);
      t300 = YY[3]*XX[4];
      t301 = XX[3]*YY[4];
      t304 = xi[0]+1.0;
      t334 = -XX[1]*t274;
      t336 = (-XX[2]-XX[3]+XX[5])*YY[1];
      t344 = -XX[1]*t281;
      t346 = YY[1]*(t27+XX[3]-XX[5]);
      t382 = t86+t262-t263+t264;
      t405 = -xi[1]+1.0;
      t448 = xi[0]-1.0;
      Pre_K_output[0][0] = t49*t4*t3;
      Pre_K_output[0][1] = t88*t4*t53;
      MapleGenVar1 = t164;
      MapleGenVar3 = 1/(t24*(-XX[1]*t55-XX[2]*t57+(XX[0]-XX[3]+XX[5])*YY[1]+(-
XX[0]+XX[3]-XX[4])*YY[2]+t66+t141+(YY[5]+YY[0])*XX[4]+(-YY[4]-YY[0])*XX[5]-t11)
+xi[2]*(-XX[1]*t70-XX[2]*t73+YY[1]*(-t26+XX[3]-XX[5])+YY[2]*(t26-XX[3]+XX[4])-
XX[4]*YY[0]+XX[5]*YY[0]+t11)-XX[1]*t81+t159+t86);
      MapleGenVar4 = t24*(XX[1]*t90+XX[2]*t92-t94+t95+t96-t97)+xi[2]*(XX[1]*(-
t69+t101+t102+t103+YY[5])+XX[2]*(t72+t106-t102-t103-YY[4])+YY[1]*(t109-t110-
t112)+YY[2]*(t115+t110+t116+XX[4])+XX[3]*t10*t100+YY[3]*t100*t16+t124*(xi[0]+xi
[1]))+XX[1]*(YY[2]+t106-t102-t103)+XX[2]*(-YY[1]+t101+t102+t103)-t85*t132;
      MapleGenVar2 = MapleGenVar3*MapleGenVar4;
      Pre_K_output[0][2] = MapleGenVar1*MapleGenVar2;
      Pre_K_output[0][3] = -t189*t167*t4;
      Pre_K_output[0][4] = -t199*t4*t192;
      Pre_K_output[0][5] = t164/(t24*(t239+t240+t241+t243+XX[3]*t170+(XX[1]-XX
[4]+XX[5])*YY[3]-XX[5]*t92-t248*YY[5])+xi[2]*(XX[1]*YY[5]+t252+t253+t254+t256+
t257-t258-t96)+t86+t262-t263+t264)*(t24*(YY[0]*t202-XX[0]*t90+t204-t205+t206-
t207)+xi[2]*(XX[0]*(-YY[5]*t111-t102+t106+t69)+YY[0]*(-t27+t115+t110+t112)+XX
[2]*(YY[3]*t215+t102+t103)+YY[2]*(-XX[3]*t215-t110-t116)+XX[3]*(xi[0]*t10+YY[5]
)+YY[3]*(xi[0]*t16-XX[5])-t124*xi[0])+XX[0]*(-YY[2]+t101+t102+t103)+YY[0]*(XX
[2]+t109-t110-t116)+XX[2]*(t106-t102-t103)-t132*YY[2]);
      Pre_K_output[0][6] = t287*t270*t4;
      Pre_K_output[0][7] = t296*t4*t290;
      Pre_K_output[0][8] = t164/(t24*(t239+t240+t334+t336+XX[3]*t272+(-XX[2]-XX
[4]+XX[5])*YY[3]-XX[4]*t90+t202*YY[4])+xi[2]*(-XX[2]*YY[4]-t204+t205+t252+t253+
t344+t346+t94)+t86+t262-t263+t264)*(t24*(-YY[0]*t248-XX[0]*t92-t257+t258-t300+
t301)+xi[2]*(XX[0]*(YY[4]*t304+t101+t103-t72)+YY[0]*(-XX[4]*t304+t109-t116+t30)
+XX[1]*(-YY[3]*t215-t102-t103)+YY[1]*(XX[3]*t215+t110+t116)+XX[3]*(xi[1]*t10-YY
[4])+YY[3]*(xi[1]*t16+XX[4])-t124*xi[1])+XX[0]*(YY[1]+t106-t102-t103)+YY[0]*(-
XX[1]+t115+t110+t116)+XX[1]*(t101+t102+t103)+t132*YY[1]);
      Pre_K_output[0][9] = -t49*xi[2]*t3;
      Pre_K_output[0][10] = -t88*xi[2]*t53;
      Pre_K_output[0][11] = t164/(t24*(-XX[4]*t62+t239+t240+t241+t292+t334+t66+
t77+t78)+xi[2]*(t252+t253+t344+t254+t60-t78)+t86+t262-t263+t264)*(t24*(-XX[1]*
t90-XX[2]*t92+t94-t95-t96+t97)+xi[2]*(-XX[0]*t1*t100-YY[0]*t100*t51+XX[1]*(-xi
[0]*t10+YY[2]*t215-YY[5])+XX[2]*(-xi[1]*t10-YY[1]*t215+YY[4])+YY[1]*(-xi[0]*t16
+XX[5])+YY[2]*(-xi[1]*t16-XX[4]))-t382*t100);
      Pre_K_output[0][12] = t189*t167*xi[2];
      Pre_K_output[0][13] = t199*xi[2]*t192;
      Pre_K_output[0][14] = t164/(t24*(t239+t240+t241+t243+t182+t77+t66+t141-
t97+t95)+xi[2]*(t252+t253+t254+t256+t172+t60)+t86+t262-t263+t264)*(t24*(-YY[0]*
t202+XX[0]*t90-t204+t205-t206+t207)+xi[2]*(XX[0]*(-YY[2]*t304+(YY[1]-YY[3]+YY
[5])*xi[0]+YY[3]*t405+t103)+YY[0]*(XX[2]*t304+(-XX[1]+XX[3]-XX[5])*xi[0]-XX[3]*
t405-t116)+XX[2]*(-YY[3]*t405-xi[0]*YY[1]-t103)+YY[2]*(XX[3]*t405+xi[0]*XX[1]+
t116)-xi[0]*(t172+t60))-xi[0]*t382);
      Pre_K_output[0][15] = -t287*xi[2]*t270;
      Pre_K_output[0][16] = -t296*xi[2]*t290;
      Pre_K_output[0][17] = t164/(t24*(t239+t240+t334+t336+t283+t78+t66+t141-
t97+t95)+xi[2]*(t252+t253+t344+t346+t276-t78)+t86+t262-t263+t264)*(t24*(YY[0]*
t248+XX[0]*t92+t257-t258+t300-t301)+xi[2]*(XX[0]*(YY[1]*t111+(-YY[2]+YY[3]-YY
[4])*xi[1]+YY[3]*t448-t102)+YY[0]*(-XX[1]*t111+(XX[2]-XX[3]+XX[4])*xi[1]-XX[3]*
t448+t110)+XX[1]*(-YY[3]*t448+xi[1]*YY[2]+t102)+YY[1]*(XX[3]*t448-xi[1]*XX[2]-
t110)+(t283+t78)*xi[1])-xi[1]*t382);

    }


	/* The gateway routine */
    void mexFunction(int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray*prhs[])

    {
        double *XX, *YY, *Z, *xi, *Pre_K_bi;

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
		plhs[0] = mxCreateDoubleMatrix(3,6,mxREAL);

		/* Assign pointers to each input and output. */
		Pre_K_bi = mxGetPr(plhs[0]);

		/* Call the C subroutine */
        fun_K_nabla_x_u(XX,YY,Z,xi,Pre_K_bi);
	}
