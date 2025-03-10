function [C] = SM_skin(E_1, E_2, E_3, G_12, G_13, G_23, nu_12, nu_13, nu_23)
nu_21 = nu_12;
nu_31 = nu_13;
nu_32 = nu_23;

% Abaqus method
gamma = 1/(1 - nu_12*nu_21 - nu_23*nu_32 - nu_31*nu_13 - 2*nu_21*nu_32*nu_13);
D1111 = E_1*(1 - nu_23*nu_32)*gamma;
D2222 = E_2*(1 - nu_13*nu_31)*gamma;
D3333 = E_3*(1 - nu_12*nu_21)*gamma;
D1122 = E_1*(nu_21 + nu_31*nu_23)*gamma;
D1133 = E_1*(nu_31 + nu_21*nu_32)*gamma;
D2233 = E_2*(nu_32 + nu_12*nu_31)*gamma;
D1212 = G_12;
D1313 = G_13;
D2323 = G_23;

D2211 = E_2*(nu_12 + nu_32*nu_13)*gamma;
D3311 = E_3*(nu_13 + nu_12*nu_23)*gamma;
D3322 = E_3*(nu_23 + nu_21*nu_13)*gamma;

DD = [D1111     D1122   D1133   0       0       0;...
      0         D2222   D2233   0       0       0;...
      0         0       D3333   0       0       0;...
      0         0       0       D1212   0       0;...
      0         0       0       0       D1313   0;...
      0         0       0       0       0       D2323];
C = DD;
end

