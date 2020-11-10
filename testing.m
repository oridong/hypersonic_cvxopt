clear all, close all, clc

% 
% n = 20;
% N = 15;
% A = eye(n);
% for j = 1:N-1
%     j
%     ind_start = (j-1)*n + 1;
%     ind_stop = j*n;
%     A_c(ind_start:ind_stop,:) = A;
% %end
% 
% [m, n] = size(A_c);
% 
% N = m/n + 1;

% [num denom] = numden(matrix)
 
% Set number of temporal nodes N if not given
addpath('matfiles/')
load('ev_test.mat')

% Set values
n = ev.opt_in.n;
N = ev.opt_in.N;
dt = ev.opt_in.dt;
x0 = ev.opt_in.x0;
u0 = ev.opt_in.u0;

% Create symbolic things
%[J_sym, dJ_sym] = ev.cost_sym();
%[A_sym, B_sym]  = ev.linsys_sym();

% Pull out linearized matrices and cost, discretize
% [A_c, B_c, fx_c] = ev.linsys_c(ev.opt_in.A_sym,ev.opt_in.B_sym,x0,u0,N,n);    
% [A_d, B_d, fx_d] = ev.linsys_d(A_c,B_c,fx_c);
% [J_c, dJ_c] = ev.cost_c(ev.opt_in.J_sym,ev.opt_in.dJ_sym,x0);
% [J_d, dJ_d] = ev.cost_d(J_c,dJ_c);

A_c = ev.opt_in.A_c;
B_c = ev.opt_in.B_c;
fx_c = ev.opt_in.fx_c;
A_d = ev.opt_in.A_d;
B_d = ev.opt_in.B_d;
fx_d = ev.opt_in.fx_d;
J_c = ev.opt_in.J_c;
dJ_c = ev.opt_in.dJ_c;
J_d = ev.opt_in.J_d;
dJ_d = ev.opt_in.dJ_d;

% Trust region constraints
    for j=1:N
        I.eps_conv(3) = abs(1/x0(3,j));
        for ind = 1:n
           norm(x(n,j) - x0(n,j)) <= I.delta_tr(n);
           norm(x(n,j) - x0(n,j))
        end
    end

