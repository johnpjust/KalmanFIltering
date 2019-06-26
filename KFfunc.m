function x = KFfunc(z,Q,Pminit, xminit, H, R, Phi, x)
    I = eye(length(Phi));
    Pm=Pminit;
    xm=xminit;

    for k=1:length(z)
        K=Pm*H'*(H*Pm*H' + R)^(-1);
        zm = H*xm;
        xhatk=xm + K*(z(k) - zm);
        x(k,:) = xhatk;
        P=(I - K*H)*Pm;
        xm=Phi*xhatk;
        Pm=Phi*P*Phi' + Q;
    end

end
