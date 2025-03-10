function [T1 T2,ilo,ihi] = balance2(A,B,rho)
%BALANCE2   Balancing generalized eigenvalue problem
%   [T1 T2] = BALANCE(A,B) computes matrices T1 and T2 that balance
%   matrices A and B for the generalized eigenvalue problem by T1*A*T2 and
%   T1*B*T2 according to [War1981].
%
%   [T1 T2] = BALANCE(A,B,rho) uses rho as computer radix. The default is
%   rho = 2.
%
%   [T1 T2 ILO IHI] = BALANCE(A,B) returns the indices ILO and IHI.
%
%
%   Example 1, see [War1981, pp.148f.]:
%
%       A = [-2.0e+1  -1.0e+4  -2.0e+0  -1.0e+6  -1.0e+1  -2.0e+5
%             6.0e-3   4.0e+0   6.0e-4   2.0e+2   3.0e-3   3.0e+1
%            -2.0e-1  -3.0e+2  -4.0e-2  -1.0e+4   0.0e+0   3.0e+3
%             6.0e-5   4.0e-2   9.0e-6   9.0e+0   3.0e-5   5.0e-1
%             6.0e-2   5.0e+1   8.0e-3  -4.0e+3   8.0e-2   0.0e+0
%             0.0e+0   1.0e+3   7.0e-1  -2.0e+5   1.3e+1  -6.0e+4 ];
%   
%       B = [-2.0e+1  -1.0e+4   2.0e+0  -2.0e+6   1.0e+1  -1.0e+5
%             5.0e-3   3.0e+0  -2.0e-4   4.0e+2  -1.0e-3   3.0e+1
%             0.0e+0  -1.0e+2  -8.0e-2   2.0e+4  -4.0e-1   0.0e+0
%             5.0e-5   3.0e-2   2.0e-6   4.0e+0   2.0e-5   1.0e-1
%             4.0e-2   3.0e+1  -1.0e-3   3.0e+3  -1.0e-2   6.0e+2
%            -1.0e+0   0.0e+0   4.0e-1  -1.0e+5   4.0e+0   2.0e+4 ];
%
%       [T1 T2] = balance2(A,B);
%       lambda1 = sort(eig(A,B)); lambda2 = sort(eig(T1*A*T2,T1*B*T2));
%       fprintf(1,'%-16s     %-16s\n','unbalanced','balanced');
%       for i = 1:size(A,1),
%           fprintf(1,'%+16.15f   %+16.15f\n',lambda1(i),lambda2(i));
%       end
%
%
%   Example 2, see [War1981, pp.148f.]:
%
%       A = full(spdiags(10.^[1:2:13; 0:2:12; -1:2:11]',[-1:1],7,7));
%       B = eye(7);
%
%       [T1 T2] = balance2(A,B);
%       lambda1 = eig(A,B); lambda2 = eig(T1*A*T2,T1*B*T2);
%       fprintf(1,'%-16s         %-16s\n','unbalanced','balanced');
%       for i = 1:size(A,1),
%           fprintf(1,'%+16.15e   %+16.15e\n',lambda1(i),lambda2(i));
%       end
%           
%
%   References:
%   [War1981]   Ward, R.C.: "Balancing the generalized eigenvalue problem",
%               SIAM Journal on Scientific and Statistical Computing, Vol.
%               2, No. 2, pp.141-152, 1981.
%
%   See also BALANCE and EIG.

%   $Version: 0.0$  $Date: 2015/02/10 08:58 $
%   Author: Mario Weder (weder@imes.mavt.ethz.ch/01-162-338)
%   Date: 2015-02-10


% Check input arguments
if (nargin < 2)
    error(message('Gcg:NotEnoughInputs'));
end

if ~isequal(size(A),size(B)),
    error('Balance2:nonsymmetric',...
        'Matrices A and B must both be square and of equal size');
end


% Default values
if nargin < 3 || isempty(rho),
    rho = 2;
end



% -------------------------------------------------------------------------
% Step 1: Reduce the order

n = size(A,1);                          % matrix size

% Incidence and permutaiton matrices
W = (A ~= 0 | B ~= 0);                  % incidence matrix, see [War1981, p.144]
p1 = 1:n;                               % incidence vector for row permuation matrix
p2 = 1:n;                               % incidence vector for column permuation matrix

% Row scanning, see [War1981, p.144]
ihi = n;
for i = n:-1:1,
    r = find(sum(W(1:i,1:i),2) == 1,1); % find row with one nonzero element
    if isempty(r),
        break;
    else
        c = find(W(r,1:i));
        p1([r i]) = p1([i r]);          % swap rows
        p2([c i]) = p2([i c]);          % swap columns
        W([r i],:) = W([i r],:);        % swap rows
        W(:,[c i]) = W(:,[i c]);        % swap columns
        ihi = i - 1;                    % update high index
    end
end

% Column scanning, see [War1981, p.144]
ilo = 1;
for i = 1:ihi,
    c = find(sum(W(i:ihi,i:ihi),1) == 1,1) + i - 1; % find column with one nonzero element
    if isempty(c),
        break;
    else
        r = find(W(i:ihi,c)) + i - 1;
        p1([r i]) = p1([i r]);          % swap rows
        p2([c i]) = p2([i c]);          % swap columns
        W([r i],:) = W([i r],:);        % swap rows
        W(:,[c i]) = W(:,[i c]);        % swap columns
        ilo = i + 1;                    % update low index
    end
end

P1 = sparse(1:n,p1,ones(1,n),n,n);      % row permutation matrix
P2 = sparse(p2,1:n,ones(1,n),n,n);      % column permuation matrix

% Return permutation matrices if no scaling is necessary
if ihi == ilo,
    T1 = P1; T2 = P2;
    ihi = 1; ilo = n;
    return;
end



% -------------------------------------------------------------------------
% Step 2: Scale A and B

% Non-triangularizable submatrices
m = ihi-ilo+1;                          % submatrix size
A22 = P1(ilo:ihi,:)*A*P2(:,ilo:ihi);    % submatrix of A
B22 = P1(ilo:ihi,:)*B*P2(:,ilo:ihi);    % submatrix of B

% Incidence matrices of nonzero elements
IA = (A22 ~= 0);                        % incidence matrix for A
IB = (B22 ~= 0);                        % incidence matrix for B

% Diagonal matrices whoes elements are the number of nonzeros, see
% [War1981, p.144]
F1 = diag(sum(IA+IB,2));                % number of nonzeros in the rows
F2 = diag(sum(IA+IB,1));                % number of nonzeros in the columns

% Sum of the incidence matrices, see [War1981, p.144]
E = IA + IB;

% Vectors of row and column sums of logarithms of nonzero elements, see
% [War1981, p.144]
logA = log(abs(A22))/log(rho); logA(~IA) = 0;
logB = log(abs(B22))/log(rho); logB(~IB) = 0;
g = sum(logA+logB,2);                   % row sums
h = sum(logA+logB,1)';                  % column sums

% Coefficient matrix to solve for powers of rho for scaling factors, see
% [War1981, p.144]
MN = [ F1  E
       E'  F2 ];

% Right hand side vector, see [War1981, p.144]
b = [ -g
      -h ];

% Inverse of auxiliary system [War1981, p.145]
I = eye(m);
e = ones(m,1);
Minv = 1/2*[1/m*I - 3/(4*m^2)*(e*e')   1/(4*m^2)*(e*e')
            1/(4*m^2)*(e*e')           1/m*I - 3/(4*m^2)*(e*e')];

% Find powers of rho by generalized conjugate gradient method using Minv
% as auxiliary matrix
rc = round(gcg(MN,b,.5,m,Minv));

% Powers of rho for row and columns scaling
r = rc(1:m);
c = rc(m+1:end);

% Global scaling matrix
D1 = speye(n); D1(ilo:ihi,ilo:ihi) = diag(rho.^r);
D2 = speye(n); D2(ilo:ihi,ilo:ihi) = diag(rho.^c);



% -------------------------------------------------------------------------
% Step 3: Grade A/B

A22 = diag(rho.^r)*A22*diag(rho.^c);    % scaled submatrix of A
B22 = diag(rho.^r)*B22*diag(rho.^c);    % scaled submatrix of B

% ratio of row and columns norms, see [War1981, p.146]
rnorm = sum(abs(A22),2)./sum(abs(B22),2);   % ratio of row norms
cnorm = sum(abs(A22),1)./sum(abs(B22),1);   % ratio of column norms

[rnorm ir] = sort(rnorm,'descend');
[cnorm ic] = sort(cnorm,'descend');

% Permutation matrices
G1 = speye(n); G1(ilo:ihi,:) = G1(ilo+ir-1,:);
G2 = speye(n); G2(:,ilo:ihi) = G2(:,ilo+ic-1);



% -------------------------------------------------------------------------
% Output

T1 = G1*D1*P1;                          % left multiplication matrix
T2 = P2*D2*G2;                          % right multiplication matrix