function [p, C] = SM_HC(t, l, h, b, theta, p_s, E_s, G_s, nu_s)
    % Density
    p = (((t/l)*((h/l)+1))*p_s)/...
        (((h/l)+sinpi(theta))*cospi(theta));

    % Elastic modulus in X direction
    E_1 = (((t/l)^3)*cospi(theta)*E_s)/...
          (((h/l) + sinpi(theta))*(sinpi(theta))^2);
    % Elastic modulus in Y direction
    E_2 = (((t/l)^3)*((h/l) + sinpi(theta))*E_s)/...
          ((cospi(theta))^3);
    % Elastic modulus in Z direction
    E_3 = ((t/l)*(1 + (h/l))*E_s)/...
          (((h/l) + sinpi(theta))*cospi(theta));
    % Shear modulus 12
    G_12 = ((t/l)^3*((h/l)+sinpi(theta))*E_s)/...
           ((h/l)^2*cospi(theta)*(1+16*(h/l)));
    % Shear modulus 13
    G_13 = ((t/l)*cospi(theta)*G_s)/...
           ((h/l)+sinpi(theta));
    % Shear modulus 31;
    G_31 = G_13;
    % Shear modulus 23
    % G_23_l = ((t/l)*((h/l)+sinpi(theta)))*G_s/...
    %          (((h/l)+1)*cospi(theta));
    % G_23_u = ((t/l)*((h/l)+(sinpi(theta))^2))*G_s/...
    %          (((h/l)+1)*cospi(theta));
    % It should be proved later, because of misunderstanding
    G_23_l = (1+(h/l)*sinpi(theta)*(t/l))*G_s/...
             ((h/l)*(2 + (h/l))*cospi(theta));
    G_23_u = ((1+2*(h/l)*(sinpi(theta))^2)*(t/l))*G_s/...
             (2*(h/l)*(1+(h/l)*sinpi(theta))*cospi(theta));
    G_23 = G_23_l + (0.787/(b/l))*(G_23_u - G_23_l);

    % Poisson's ratios
    % nu_12 = ((cospi(theta))^2)/...
    %         (((h/l) + sinpi(theta))*sinpi(theta));
    nu_12 = 0.99;
    % nu_21 = (((h/l) + sinpi(theta))*sinpi(theta))/...
    %         ((cospi(theta))^2);
    nu_21 = 0.99;
    nu_31 = nu_s;
    nu_32 = nu_s;
    nu_13 = nu_s*(E_1/E_3);
    nu_23 = nu_s*(E_2/E_3);

    % Stiffness matrix
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

