%%% Created by Amirpasha Shirazinia (November 4, 2014), Email:
%%% amirpasha.shirazinia@signal.uu.se
%%% You need to instal CVX first, and then ctrl(or command)+enter each
%%% cell.

%% Scalar example -- eq. (2)

clc
clear all

N = 10; alpha = rand(1,N);
cvx_begin
    variable p(N) 
    p >= 0
    sum(p) == 1
    minimize(alpha*p)
cvx_end

%% Convexifying a non-convex problem -- eq. (5)
clc
clear all
N = 10;
cvx_begin 
    variables x(N,1) y z
    y >= 0; z>= 0;
    norm([2*x ; y-z],2) <= y + z
    minimize(x'*x + y + z)
cvx_end

%% Measurement matrix design in a Bayesian (MMSE) framework -- eq. (9)
clc
clear all

N = 3; rho = 0.5; P = 10; sigma_sq = 0.01; 
R = [1 rho rho^2; rho 1 rho; rho^2 rho 1];
cvx_begin sdp
    variable X(N,N) hermitian 
    X == semidefinite(N,N)
    %X >= 0
    trace(R*X) <= P
    minimize(trace_inv(R^-1 + X/sigma_sq))
cvx_end
A = chol(X);

%% l1-norm relaxation in compressed sensing framework -- eq. (11)
clc
clear all
close all

N = 50; M = 25; K = 8;
supp = randsample(N,K)'; x = zeros(N,1); x(supp) = randn(K,1);
A = randn(M,N); s = sqrt(sum(A.^2)); S = diag(1./s);
A = A*S;
y = A*x;
cvx_begin
    variable x_hat(N,1)
    y == A*x_hat
    minimize (norm(x_hat,1))
cvx_end
[x x_hat]
plot(1:N,x,'*',1:N,x_hat,'s')
