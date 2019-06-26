% PROGRAM NAME: linekfwn.m
% PURPOSE: Track a t.v. line with assumed WN contamination
% REQUIRED INPUT: from program linegenwn.m include:
% z = time series including line + WN
% npts=length(z)
% sn2 = additive WN variance
% =======================================================
% npts = length(SCALE_1)-1;
% z = SCALE_1(2:end);
dims = 3;
npts = length(z);
I = eye(dims);
Phi= [1 0 0; 1 1 0;0 0 ar1];
% Cov Matrix for x is obtained from trial+ error:
sa2= 10.^-(1:1:15); % <==== SLOPE DRIVING NOISE VARIANCE =====to be entered
% sa2= [10.^1 10.^-5];
% The following is the additive WN variance from linegenwn.m:
R=0;
%=====================================
% EKF Initial Conditions
%====================================


x = zeros(length(sa2),npts,dims);
mse = zeros(length(sa2),1);
% EKF Loop
for nvar = 1:length(sa2)
    Q=[sa2(nvar) 0 0;0 0 0;0 0 ns2];
    Pm=Q;
    xm=[0; 0; 0]; 
    temp_shat = zeros(dims,npts);
    H=[0 1 1];
    for k=1:npts
        
        K=Pm*H'*(H*Pm*H' + R)^(-1);
        zm = H*xm;
        xhatk=xm + K*(z(k) - zm);
        x(nvar,k,:) = xhatk;
        temp_shat(:,k) = xhatk;
        P=(I - K*H)*Pm;
        xm=Phi*xhatk;
        Pm=Phi*P*Phi' + Q;
    end
%     mse(nvar) = mean((z-H*temp_shat).^2);
    mse(nvar) = mean((s-[0 1 0]*temp_shat).^2);
end
% ind = nvar;
[val, ind] = min(mse);
    ahat = x(ind,:,1);
    shat=x(ind,:,2);
% =================================================================

%=====================
tvec=1:npts;
figure(1)
plot(tvec,z,tvec,s,tvec,shat,'k')
title('Comparison of Raw Data and KF-Estimate of T.V. Line')
% 
% plot(tvec,z,'b',tvec,x(ind,:,2),'g', tvec,x(ind,:,2),'r')
% title('KF-Estimate of T.V. Line')
grid

figure(2)
clf
plot(tvec,av,'LineWidth',2);hold on; plot(tvec,x(ind,:,1),'LineWidth',2)
% plot(tvec,x(ind,:,1),'Color','r','LineWidth',2)
title('KF Estimate of T.V. Line Slope')
xlabel('Time [sec]')
ylabel('Slope')
grid

% figure(3) 
% plot(tvec,n,tvec,x(ind,:,3))

% pause
% end
figure(3);
loglog(sa2,mse)
