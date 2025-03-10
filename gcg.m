function [x k] = gcg(A,c,tol,maxit,Minv,x0)
%GCG   Generalized conjugate gradient method
%   X = GCG(A,B)  attempts to solve the system of linear equations A*X=B
%   for X. The N-by-N coefficient matrix A must be symmetric and positive
%   definite and the right hand side column vector B must have length N.
%  
%   X = GCG(A,B,TOL) specifies the tolerance of the method. If TOL is []
%   then GCG uses the default, 1e-6.
%
%   X = GCG(A,B,TOL,MAXIT) specifies the maximum number of iterations. If
%   MAXIT is [] then GCG uses the default, min(N,20).
%
%   X = GCG(A,B,TOL,MAXIT,MINV) provides auxiliary matrix inverse Minv of a
%   simpler solvable system. If MINV is [] then GCG uses the default,
%   eye(N,N).
%
%   X = GCG(A,B,TOL,MAXIT,M,X0) specifies the initial guess. If X0 is []
%   then GCG uses the default, an all zero vector.
%
%   [X K] = GCG(A,B) returns the number of iterations K.
%
%
%   References:
%   [CGO1976]   Concus, P., Golub, G.H. and O'Leary D.P.: "A generalized
%               conjugate gradient method for the numerical solution of
%               elliptic partial differential equations", in Sparse Matrix
%               Computations, Bunch J. R and Rose, D. J., eds., Academic
%               Press, 1976, pp. 309-332.
%
%   See also PCG, BICG. BICGSTAB, BICGSTABL and CGS.

%   $Version: 0.0$  $Date: 2015/02/10 08:58 $
%   Author: Mario Weder (weder@imes.mavt.ethz.ch/01-162-338)
%   Date: 2015-02-10


% Check input
if (nargin < 2)
    error(message('Gcg:NotEnoughInputs'));
end

% Check matrix and right hand side vector inputs have appropriate sizes
[n,m] = size(A);
if (n ~= m)
    error(message('Gcg:NonSquareMatrix'));
end
if ~isequal(size(c),[n,1])
    error(message('Gcg:RSHsizeMatchCoeffMatrix',n));
end

    
% Default values
if nargin < 3 || isempty(tol),
    tol = 1e-6;
end

if nargin < 4 || isempty(maxit),
    maxit = min(n,20);
end

if nargin < 5 || isempty(Minv),
    Minv = eye(n,n);
end

if nargin < 6 || isempty(x0),
    x0 = zeros(n,1);
end


% Initial guess
p = x0;
x = x0;

for k = 1:maxit,
    
    % Step 1, see [CGO1975, p.315]
    Mz = c - A*x;
    z = Minv*Mz;
    
    % Step 2, see [CGO1975, p.315]
    if k == 1,
        b = 0;
    else
        b = z'*Mz/zMz;                  % use zMz of previous step
    end
    p = z + b*p;
    
    % Step 3, see [CGO1975, p.316]
    zMz = z'*Mz;                        % update zMz
    a = zMz/(p'*A*p);
    x = x + a*p;
    
    % Break if update increment is below tolerance
    if norm(a*p) < tol, break; end
end