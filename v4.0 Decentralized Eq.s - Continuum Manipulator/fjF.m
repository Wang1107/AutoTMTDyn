function fj = fjF(lambdal,in2)
%FJF
%    FJ = FJF(LAMBDAL,IN2)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    24-Sep-2017 19:43:37

q1 = in2(:,1);
q2 = in2(:,2);
q3 = in2(:,3);
u1 = in2(:,4);
u2 = in2(:,5);
u3 = in2(:,6);
fj = [u1.*(-1.0./1.0e3)-1.0./lambdal.^2.*(q1-1.1e1./2.5e2).*1.406576710811794e3;u2.*(-4.0e-5)-1.0./lambdal.^3.*(q2-1.0./1.0e3).*2.031373444765246e-3;u3.*(-4.0e-5)-1.0./lambdal.^3.*(q3-1.0./1.0e3).*2.031373444765246e-3];
