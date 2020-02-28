function [Hinf,winf3] = Inf_norm_discrete(A,B,C,D,wl,wu,tol,ts)
Ac = -(inv((eye(length(A))+A))) * (eye(length(A))-A);
Bc = sqrt(2) * (inv((eye(length(A))+A))) * B;
Cc = sqrt(2) * C * (inv((eye(length(A))+A)));
Dc = D - C * (inv((eye(length(A))+A))) * B;
[Hinf,winf1] = Inf_norm_continuous(Ac,Bc,Cc,Dc,wl,wu,tol);
winf2 = (log((1+(1i*winf1))/(1-(1i*winf1))))/(1i*ts);
winf3 = real(winf2);
end

