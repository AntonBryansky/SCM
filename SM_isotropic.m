function [C] = SM_isotropic(E, nu)

mu = E/(2*nu);
lambda = E*nu/((1 + nu)*(1 - 2*nu));

% Abaqus method
D1111 = lambda + 2*mu;
D1122 = lambda;
D1133 = lambda;
D2222 = lambda + 2*mu;
D2233 = lambda;
D3333 = lambda + 2*mu;
D1212 = mu;
D1313 = mu;
D2323 = mu;

C = [D1111     D1122   D1133   0       0       0;...
     0         D2222   D2233   0       0       0;...
     0         0       D3333   0       0       0;...
     0         0       0       D1212   0       0;...
     0         0       0       0       D1313   0;...
     0         0       0       0       0       D2323];
end

