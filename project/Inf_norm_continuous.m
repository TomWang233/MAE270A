function [Hinf,winf] = Inf_norm_continuous(A,B,C,D,gamma_l,gamma_u,tol)
    while(gamma_u - gamma_l) > tol
    gamma = (gamma_l + gamma_u)/2;
    D_gamma = (gamma^2)*eye(size(D'*D,1)) - D'*D;
    A1 = A+B*inv(D_gamma)*D'*C;
    D1 = D*inv(D_gamma)*D';
    A_clp = [A1,(B*inv(D_gamma)*B');-C'*(eye(size(D1))+D1)*C, -A1'];
    eig_test = eig(A_clp);
    if (min(abs(real(eig_test))) < 1e-10) 
        gamma_l = gamma;
    else
        gamma_u = gamma;
    end
    Hinf = gamma;
    winf = max(imag(eig_test));
    end  
end

