function TC = TCF(lambdal,in2,in3)
%TCF
%    TC = TCF(LAMBDAL,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    29-Sep-2017 11:21:43

q1 = in2(:,1);
q2 = in2(:,2);
q3 = in2(:,3);
q4 = in2(:,4);
q6 = in2(:,6);
s1 = in3(:,1);
s2 = in3(:,2);
s3 = in3(:,3);
u2 = in2(:,8);
u4 = in2(:,10);
u6 = in2(:,12);
t2 = 1.0./q2;
t3 = cos(q2);
t4 = 1.0./q2.^2;
t5 = sin(q2);
t6 = t3-1.0;
t7 = 1.0./q4;
t8 = sin(q4);
t9 = s2.*t3.*t7.*t8;
t10 = cos(q4);
t11 = 1.0./q4.^2;
t12 = t10-1.0;
t13 = q1.*t2.*t5;
t14 = q1.*t4.*t6;
t15 = 1.0./q6;
t16 = sin(q6);
t17 = t3.*t10;
t27 = t5.*t8;
t18 = t17-t27;
t19 = s3.*t15.*t16.*t18;
t20 = q3.*t3.*t7.*t8;
t21 = cos(q6);
t22 = t21-1.0;
t23 = t3.*t8;
t24 = t5.*t10;
t25 = t23+t24;
t26 = s3.*t15.*t22.*t25;
t28 = 1.0./q6.^2;
t29 = t2.*t5;
t30 = q1.*t2.*t3;
t31 = s3.*t15.*t18.*t22;
t32 = t3.^2;
t33 = t5.^2;
t34 = t32.*u2;
t35 = t33.*u2;
t36 = t34+t35;
t37 = t25.*u2;
t38 = t25.*u4;
t39 = t37+t38;
t40 = t25.*t39;
t41 = t18.*u2;
t42 = t18.*u4;
t43 = t41+t42;
t44 = t18.*t43;
t45 = t40+t44;
t46 = t25.^2;
t47 = t18.^2;
t48 = t46+t47;
t49 = t21.*t25;
t50 = t16.*t18;
t51 = t49+t50;
t52 = t18.*t21;
t54 = t16.*t25;
t53 = t52-t54;
t55 = t51.*u2;
t56 = t51.*u4;
t57 = t51.*u6;
t58 = t55+t56+t57;
t59 = t51.*t58;
t60 = t53.*u2;
t61 = t53.*u4;
t62 = t53.*u6;
t63 = t60+t61+t62;
t64 = t53.*t63;
t65 = t59+t64;
t66 = t51.^2;
t67 = t53.^2;
t68 = t66+t67;
TC = reshape([0.0,0.0,0.0,-t2.*t6,0.0,t29,-t2.*t6,0.0,t29,0.0,t36,0.0,0.0,t45,0.0,0.0,t65,0.0,s1.*t2.*t5+s1.*t4.*t6,0.0,s1.*t2.*t3-s1.*t4.*t5,t9+t13+t14+s2.*t5.*t7.*t12,0.0,t30-q1.*t4.*t5-s2.*t5.*t7.*t8+s2.*t3.*t7.*t12,t13+t14+t19+t20+t26+q3.*t5.*t7.*t12,0.0,t30+t31-q1.*t4.*t5-q3.*t5.*t7.*t8+q3.*t3.*t7.*t12-s3.*t15.*t16.*t25,0.0,t32+t33,0.0,0.0,t48,0.0,0.0,t68,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t7.*t8-t3.*t7.*t12,0.0,t3.*t7.*t8+t5.*t7.*t12,0.0,t36,0.0,0.0,t45,0.0,0.0,t65,0.0,0.0,0.0,0.0,t9+s2.*t5.*t7.*t10-s2.*t5.*t8.*t11+s2.*t3.*t11.*t12,0.0,s2.*t3.*t7.*t10-s2.*t5.*t7.*t8-s2.*t3.*t8.*t11-s2.*t5.*t11.*t12,t19+t20+t26+q3.*t5.*t7.*t10-q3.*t5.*t8.*t11+q3.*t3.*t11.*t12,0.0,t31+q3.*t3.*t7.*t10-q3.*t5.*t7.*t8-q3.*t3.*t8.*t11-q3.*t5.*t11.*t12-s3.*t15.*t16.*t25,0.0,t36,0.0,0.0,t48,0.0,0.0,t68,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t36,0.0,0.0,t45,0.0,0.0,t65,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t19+s3.*t15.*t21.*t25+s3.*t18.*t22.*t28-s3.*t16.*t25.*t28,0.0,s3.*t15.*t18.*t21-s3.*t15.*t16.*t25-s3.*t16.*t18.*t28-s3.*t22.*t25.*t28,0.0,t36,0.0,0.0,t45,0.0,0.0,t68,0.0],[18,6]);