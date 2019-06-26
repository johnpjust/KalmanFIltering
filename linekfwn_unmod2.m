% PROGRAM NAME: linekfwn.m
% PURPOSE: Track a t.v. line with assumed WN contamination
% REQUIRED INPUT: from program linegenwn.m include:
% z = time series including line + WN
% npts=length(z)
% sn2 = additive WN variance
% =======================================================
I = eye(2,2);
Phi= [1 0  ; 1 1];
% Cov Matrix for x is obtained from trial+ error:
sa2= 10.^-(12) % <==== SLOPE DRIVING NOISE VARIANCE =====to be entered
Q=[sa2 0 ; 0 0];
% The following is the additive WN variance from linegenwn.m:
R=ns2;
%=====================================
% EKF Initial Conditions
xm=[0; 0]; 
Pm=Q;
%====================================
x = [];
% EKF Loop
for k=1:npts
H=[0 1];
K=Pm*H'*(H*Pm*H' + R)^(-1);
zm = H*xm;
xhatk=xm + K*(z(k) - zm);
x = [x xhatk];
P=(I - K*H)*Pm;
xm=Phi*xhatk;
Pm=Phi*P*Phi' + Q;
end
% =================================================================
ahat=x(1,:);
shat=x(2,:);
%=====================
tvec=1:npts;
figure(1)
plot(tvec,z,tvec,s,tvec,shat,'k')
title('Comparison of Raw Data and KF-Estimate of T.V. Line')
grid
mse = mean((s-shat).^2);
% PLOT LINE SLOPE AND INTERCEPT ESTIMATES
figure(2)
plot(tvec,av,'LineWidth',2)
hold on
plot(tvec,ahat)
title('KF Estimate of T.V. Line Slope')
xlabel('Time [sec]')
ylabel('Slope')
grid
