 	#include "mex.h"
	#include "math.h"

	void rectangule_K_gen(double X[], double Y[], double Z[], double lambda[], double K[][4])
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
    double t500, t501, t502, t503, t504, t505, t506, t507, t508, t509;
    double t510, t511, t512, t513, t514, t515, t516, t517, t518, t519;
    double t520, t521, t522, t523, t524, t525, t526, t527, t528, t529;
    double t530, t531, t532, t533, t534, t535, t536, t537, t538, t539;
	  
      t3 = X[0]*X[1];
      t4 = lambda[1]*lambda[1];
      t5 = t4*t4;
      t6 = t5*t3;
      t7 = 6.0*t6;
      t8 = X[1]*X[1];
      t9 = lambda[0]*lambda[0];
      t10 = t9*t9;
      t12 = 3.0*t10*t8;
      t13 = t9*lambda[0];
      t15 = lambda[1]*t13*t8;
      t16 = 12.0*t15;
      t17 = t9*t8;
      t18 = t4*t17;
      t19 = 18.0*t18;
      t20 = lambda[0]*t8;
      t21 = t4*lambda[1];
      t22 = t21*t20;
      t23 = 12.0*t22;
      t25 = 3.0*t5*t8;
      t26 = Y[0]*Y[1];
      t27 = t5*t26;
      t28 = 6.0*t27;
      t29 = Y[1]*Y[1];
      t31 = 3.0*t10*t29;
      t33 = lambda[1]*t13*t29;
      t34 = 12.0*t33;
      t35 = t9*t29;
      t36 = t4*t35;
      t37 = 18.0*t36;
      t38 = lambda[0]*t29;
      t39 = t21*t38;
      t40 = 12.0*t39;
      t42 = 3.0*t5*t29;
      t43 = Z[0]*Z[0];
      t44 = t9*t43;
      t45 = 3.0*t44;
      t48 = 12.0*lambda[1]*lambda[0]*t43;
      t49 = t4*t43;
      t50 = 9.0*t49;
      t51 = Z[1]*Z[1];
      t52 = t9*t51;
      t53 = 3.0*t52;
      t54 = t4*t51;
      t55 = 9.0*t54;
      t56 = t7+t12-t16+t19-t23+t25+t28+t31-t34+t37-t40+t42+t45-t48+t50+t53+t55;
      t57 = Z[0]*Z[1];
      t58 = t9*t57;
      t59 = 6.0*t58;
      t60 = t4*t57;
      t61 = 18.0*t60;
      t64 = 12.0*lambda[1]*lambda[0]*t51;
      t65 = log(lambda[0]);
      t66 = t65*t5;
      t67 = X[0]*X[0];
      t69 = 2.0*t67*t66;
      t71 = 2.0*t8*t66;
      t72 = Y[0]*Y[0];
      t74 = 2.0*t72*t66;
      t76 = 2.0*t29*t66;
      t77 = log(lambda[1]);
      t78 = t77*t5;
      t80 = 2.0*t67*t78;
      t82 = 2.0*t8*t78;
      t84 = 2.0*t72*t78;
      t86 = 2.0*t29*t78;
      t87 = t65*t4;
      t89 = 6.0*t43*t87;
      t91 = 6.0*t51*t87;
      t92 = t77*t4;
      t94 = 6.0*t43*t92;
      t96 = 6.0*t51*t92;
      t97 = t9*t3;
      t98 = t97*t87;
      t99 = 2.0*t98;
      t100 = t65*t21;
      t101 = lambda[0]*t3;
      t102 = t101*t100;
      t103 = 4.0*t102;
      t104 = -t59-t61-t64+t69+t71+t74+t76-t80-t82-t84-t86+t89+t91-t94-t96+t99-
t103;
      t106 = t9*t26;
      t107 = t106*t87;
      t108 = 2.0*t107;
      t109 = lambda[0]*t26;
      t110 = t109*t100;
      t111 = 4.0*t110;
      t112 = t97*t92;
      t113 = 2.0*t112;
      t114 = t77*t21;
      t115 = t101*t114;
      t116 = 4.0*t115;
      t117 = t106*t92;
      t118 = 2.0*t117;
      t119 = t109*t114;
      t120 = 4.0*t119;
      t121 = lambda[1]*t13;
      t122 = t121*t3;
      t123 = 6.0*t122;
      t124 = t4*t9;
      t126 = 18.0*t124*t3;
      t127 = t21*lambda[0];
      t128 = t127*t3;
      t130 = t121*t26;
      t131 = 6.0*t130;
      t133 = 18.0*t124*t26;
      t134 = t127*t26;
      t136 = lambda[0]*lambda[1];
      t138 = 24.0*t136*t57;
      t139 = t9*t67;
      t140 = t139*t87;
      t141 = 2.0*t140;
      t142 = lambda[0]*t67;
      t143 = t142*t100;
      t144 = 4.0*t143;
      t146 = 2.0*t3*t66;
      t147 = t17*t87;
      t148 = 2.0*t147;
      t149 = t108-t111-t113+t116-t118+t120-t123+t126-18.0*t128-t131+t133-18.0*
t134+t138+t141-t144+t146+t148;
      t150 = t20*t100;
      t151 = 4.0*t150;
      t152 = t9*t72;
      t153 = t152*t87;
      t154 = 2.0*t153;
      t155 = lambda[0]*t72;
      t156 = t155*t100;
      t157 = 4.0*t156;
      t159 = 2.0*t26*t66;
      t160 = t35*t87;
      t161 = 2.0*t160;
      t162 = t38*t100;
      t163 = 4.0*t162;
      t164 = t139*t92;
      t165 = 2.0*t164;
      t166 = t142*t114;
      t167 = 4.0*t166;
      t169 = 2.0*t3*t78;
      t170 = t17*t92;
      t171 = 2.0*t170;
      t172 = t20*t114;
      t173 = 4.0*t172;
      t174 = t152*t92;
      t175 = 2.0*t174;
      t176 = t155*t114;
      t177 = 4.0*t176;
      t179 = 2.0*t26*t78;
      t180 = t35*t92;
      t181 = 2.0*t180;
      t182 = t38*t114;
      t183 = 4.0*t182;
      t185 = 12.0*t57*t87;
      t187 = 12.0*t57*t92;
      t188 = -t151+t154-t157+t159+t161-t163-t165+t167-t169-t171+t173-t175+t177-
t179-t181+t183-t185+t187;
      t189 = t149+t188;
      t191 = t29*t67;
      t206 = t72*t8;
      t218 = 4.0*Y[1]*lambda[0]*lambda[1]*Y[0]*t3-2.0*t4*t26*t3-2.0*Z[1]*Z[0]*
t67-2.0*t106*t3-2.0*t136*t191-2.0*t136*t206+t4*t191+t9*t191+t4*t206+t9*t206-2.0
*t43*t3+t43*t67+t51*t67;
      t244 = -2.0*Z[1]*Z[0]*t29-2.0*Z[1]*Z[0]*t72-2.0*Z[1]*Z[0]*t8-2.0*t43*t26
-2.0*t51*t26+4.0*t57*t26+t43*t29+t51*t29-2.0*t51*t3+4.0*t57*t3+t43*t72+t43*t8+
t51*t72+t51*t8;
      t246 = sqrt(t218+t244);
      t247 = 1/t246;
      t249 = lambda[0]-lambda[1];
      t250 = t249*t249;
      t252 = 1/t250/t249;
      t256 = 3.0*t5*t67;
      t257 = 3.0*t6;
      t258 = 3.0*t15;
      t259 = 9.0*t18;
      t260 = 9.0*t22;
      t262 = 3.0*t5*t72;
      t263 = 3.0*t27;
      t264 = 3.0*t33;
      t265 = 9.0*t36;
      t266 = 9.0*t39;
      t267 = t256+t257-t258+t259-t260+t25+t262+t263-t264+t265-t266+t42+t45-t48+
t50-t59-t61+t53+t55;
      t269 = lambda[1]*t13*t67;
      t270 = 3.0*t269;
      t271 = t4*t139;
      t272 = 9.0*t271;
      t273 = t21*t142;
      t274 = 9.0*t273;
      t275 = t10*t3;
      t276 = 3.0*t275;
      t278 = lambda[1]*t13*t72;
      t279 = 3.0*t278;
      t280 = t4*t152;
      t281 = 9.0*t280;
      t282 = -t64+t69+t71+t74+t76-t80-t82-t84-t86+t89+t91-t94-t96-t270+t272-
t274+t276-t279+t281;
      t284 = t21*t155;
      t285 = 9.0*t284;
      t286 = t10*t26;
      t287 = 3.0*t286;
      t288 = 12.0*t122;
      t289 = 12.0*t128;
      t290 = 12.0*t130;
      t291 = 12.0*t134;
      t292 = -t285+t287+t99-t103+t108-t111-t113+t116-t118+t120-t288+t126-t289-
t290+t133-t291+t138+t141-t144;
      t293 = t146+t148-t151+t154-t157+t159+t161-t163-t165+t167-t169-t171+t173-
t175+t177-t179-t181+t183-t185+t187;
      t298 = t252*t247*(t267+t282+t292+t293)/6.0;
      t299 = 3.0*t49;
      t300 = 3.0*t54;
      t301 = 6.0*t60;
      t302 = 4.0*t98;
      t303 = 2.0*t102;
      t304 = 4.0*t107;
      t305 = 2.0*t110;
      t306 = -t45+t299-t53+t300+t257+t263+t59-t301-t276-t287-t302+t303-t304+
t305;
      t307 = 4.0*t112;
      t308 = 2.0*t115;
      t309 = 4.0*t117;
      t310 = 2.0*t119;
      t311 = t65*t13;
      t313 = X[1]*lambda[1]*X[0];
      t314 = t313*t311;
      t315 = 2.0*t314;
      t317 = Y[1]*lambda[1]*Y[0];
      t318 = t317*t311;
      t319 = 2.0*t318;
      t320 = t77*t13;
      t321 = t313*t320;
      t322 = 2.0*t321;
      t323 = t317*t320;
      t324 = 2.0*t323;
      t325 = t65*lambda[0];
      t327 = Z[0]*Z[1]*lambda[1];
      t329 = 12.0*t327*t325;
      t330 = t77*lambda[0];
      t332 = 12.0*t327*t330;
      t333 = 6.0*t128;
      t334 = 6.0*t134;
      t335 = t307-t308+t309-t310+t315+t319-t322-t324-t329+t332+t123-t333+t131-
t334;
      t337 = 4.0*t140;
      t338 = 2.0*t143;
      t339 = 4.0*t147;
      t340 = 2.0*t150;
      t341 = 4.0*t153;
      t342 = 2.0*t156;
      t343 = 4.0*t160;
      t344 = 2.0*t162;
      t345 = 4.0*t164;
      t346 = 2.0*t166;
      t347 = 4.0*t170;
      t348 = 2.0*t172;
      t349 = 4.0*t174;
      t350 = 2.0*t176;
      t351 = -t337+t338-t339+t340-t341+t342-t343+t344+t345-t346+t347-t348+t349-
t350;
      t352 = 4.0*t180;
      t353 = 2.0*t182;
      t354 = lambda[1]*t67;
      t355 = t354*t311;
      t356 = 2.0*t355;
      t357 = lambda[1]*t8;
      t358 = t357*t311;
      t359 = 2.0*t358;
      t360 = lambda[1]*t72;
      t361 = t360*t311;
      t362 = 2.0*t361;
      t363 = lambda[1]*t29;
      t364 = t363*t311;
      t365 = 2.0*t364;
      t366 = t354*t320;
      t367 = 2.0*t366;
      t368 = t357*t320;
      t369 = 2.0*t368;
      t370 = t360*t320;
      t371 = 2.0*t370;
      t372 = t363*t320;
      t373 = 2.0*t372;
      t374 = t43*lambda[1];
      t376 = 6.0*t374*t325;
      t377 = t51*lambda[1];
      t379 = 6.0*t377*t325;
      t381 = 6.0*t374*t330;
      t383 = 6.0*t377*t330;
      t384 = t352-t353+t356+t359+t362+t365-t367-t369-t371-t373+t376+t379-t381-
t383;
      t389 = t252*t247*(t306+t335+t351+t384)/6.0;
      t390 = 9.0*t15;
      t391 = 3.0*t22;
      t392 = 9.0*t33;
      t393 = 3.0*t39;
      t394 = t256-t12+t390-t259+t391+t262-t31+t392-t265+t393-t45+t299+t59-t301-
t53+t300;
      t395 = -t270+t272-t274-t279+t281-t285-t302+t303-t304+t305+t307-t308+t309-
t310+t315+t319;
      t397 = -t322-t324-t329+t332-t337+t338-t339+t340-t341+t342-t343+t344+t345-
t346+t347-t348;
      t398 = t349-t350+t352-t353+t356+t359+t362+t365-t367-t369-t371-t373+t376+
t379-t381-t383;
      t399 = t397+t398;
      t403 = t252*t247*(t394+t395+t399)/6.0;
      t405 = 3.0*t10*t67;
      t407 = 3.0*t10*t72;
      t408 = t45+t50+t53+t55+t256+t262+t405+t407+t7+t28-t48-t59-t61-t64+t69+t71
+t74;
      t409 = 12.0*t269;
      t410 = 18.0*t271;
      t411 = 12.0*t273;
      t412 = 12.0*t278;
      t413 = 18.0*t280;
      t414 = 12.0*t284;
      t415 = t76-t80-t82-t84-t86+t89+t91-t94-t96-t409+t410-t411-t412+t413-t414+
t99-t103;
      t421 = -t405-t258+t259-t260+t25-t407-t264+t265-t266+t42-t45+t299+t59-t301
-t53+t300;
      t422 = 9.0*t269;
      t423 = 3.0*t273;
      t424 = 9.0*t278;
      t425 = 3.0*t284;
      t426 = t422-t272+t423+t424-t281+t425-t302+t303-t304+t305+t307-t308+t309-
t310+t315+t319;
      t431 = t252*t247*(t421+t426+t399)/6.0;
      t432 = 9.0*t44;
      t433 = 18.0*t58;
      t434 = 9.0*t52;
      t435 = t12-t16+t19-t23+t25+t31-t34+t37-t40+t42+t432-t48+t299-t433-t301+
t434+t300;
      t436 = 6.0*t275;
      t437 = 6.0*t286;
      t438 = t65*t10;
      t440 = 2.0*t67*t438;
      t442 = 2.0*t8*t438;
      t444 = 2.0*t72*t438;
      t446 = 2.0*t29*t438;
      t447 = t77*t10;
      t449 = 2.0*t67*t447;
      t451 = 2.0*t8*t447;
      t453 = 2.0*t72*t447;
      t455 = 2.0*t29*t447;
      t456 = t65*t9;
      t458 = 6.0*t43*t456;
      t460 = 6.0*t51*t456;
      t461 = t77*t9;
      t463 = 6.0*t43*t461;
      t465 = 6.0*t51*t461;
      t466 = -t64+t436+t437-t440-t442-t444-t446+t449+t451+t453+t455-t458-t460+
t463+t465-t99-t108;
      t468 = 4.0*t314;
      t469 = 4.0*t318;
      t470 = 4.0*t321;
      t471 = 4.0*t323;
      t474 = t113+t118+t468+t469-t470-t471-18.0*t122+t126-t333-18.0*t130+t133-
t334+t138-t141-t148-t154-t161;
      t476 = 12.0*t57*t456;
      t478 = 12.0*t57*t461;
      t479 = 4.0*t355;
      t480 = 4.0*t358;
      t481 = 4.0*t361;
      t482 = 4.0*t364;
      t483 = 4.0*t366;
      t484 = 4.0*t368;
      t485 = 4.0*t370;
      t486 = 4.0*t372;
      t488 = 2.0*t3*t438;
      t490 = 2.0*t26*t438;
      t492 = 2.0*t3*t447;
      t494 = 2.0*t26*t447;
      t495 = t165+t171+t175+t181+t476-t478+t479+t480+t481+t482-t483-t484-t485-
t486-t488-t490+t492+t494;
      t496 = t474+t495;
      t501 = -t405-t257-t12+t390-t259+t391-t407-t263-t31+t392-t265+t393-t432+
t48-t299+t433+t301-t434-t300;
      t502 = t64+t422-t272+t423-t276+t424-t281+t425-t287+t440+t442+t444+t446-
t449-t451-t453-t455+t458+t460;
      t504 = -t463-t465+t99+t108-t113-t118-t468-t469+t470+t471+t288-t126+t289+
t290-t133+t291-t138+t141+t148;
      t505 = t154+t161-t165-t171-t175-t181-t476+t478-t479-t480-t481-t482+t483+
t484+t485+t486+t488+t490-t492-t494;
      t510 = t252*t247*(t501+t502+t504+t505)/6.0;
      t511 = -t405+t409-t410+t411-t256-t436-t407+t412-t262-t432+t48-t299+t433+
t301-t434+t64-t300;
      t512 = -t413+t414-t437+t440+t442+t444+t446-t449-t451-t453-t455+t458+t460-
t463-t465+t99+t108;
      K[0][0] = t252*t247*(t56+t104+t189)/6.0;
      K[0][1] = -t298;
      K[0][2] = -t389;
      K[0][3] = t403;
      K[1][0] = -t298;
      K[1][1] = t252*t247*(t408+t415+t189)/6.0;
      K[1][2] = t431;
      K[1][3] = -t389;
      K[2][0] = -t389;
      K[2][1] = t431;
      K[2][2] = -t252*t247*(t435+t466+t496)/6.0;
      K[2][3] = -t510;
      K[3][0] = t403;
      K[3][1] = -t389;
      K[3][2] = -t510;
      K[3][3] = t252*t247*(t511+t512-t496)/6.0;


    }


	/* The gateway routine */
    void mexFunction(int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray*prhs[])

    {
        double *X, *Y, *Z, *lambda, *K;
   
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
		K = mxGetPr(plhs[0]);

		/* Call the C subroutine */
        rectangule_K_gen(X,Y,Z,lambda,K);
	}
