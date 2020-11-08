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



% [c,z0,M,F] = ev.restack(x0,u0,dJ_d,A_d,B_d,fx_d,A_c,dt,N,n);


c = [transpose(dJ_d);...
                 zeros(N,1)];

x_z = reshape(x0,[n*N,1]);
u_z = reshape(u0,[N,1]);

z0 = [x_z;...
      u_z];


for j=1:N
    r1 = (j-1)*n + 1;
    r2 = j*n; 

    Adiag{1,j} = A_d(r1:r2,:);
    Bdiag{1,j} = B_d(r1:r2,:);

    Idiag{1,j} = eye(n-1);

    F(r1:r2,1) = fx_d(r1:r2,1) - dt*A_c(r1:r2,:)*x0(:,j) - B_d(r1:r2,1)*u0(j);
end

F = [x0(:,1);...
     F];
 
Ishift1 = [ eye(n*N); ...
            zeros(n,n*N)];
Zshift1 = zeros(n,n*N);
Zshift2  = zeros(n,N);

M1 = Ishift1 + [Zshift1; -blkdiag(Adiag{:})];
M2 = [  Zshift2;...
        -blkdiag(Bdiag{:})];

        
M = [M1, M2];


blah = M*z0 - F;

%size([eye(n*(N-1)),Ishift; zeros(n,N*n)])

 %{
% Reorganize cost and constraints into 
%function [c,z0,M,F] = restack(x0,u0,dJ_d,A_d,B_d,fx_d,A_c,dt,N,n)
    %if nargin==6
    %   [n, N] = size(x0);
    %   dt = opt_in.dt;
    %end 

    c = [transpose(dJ_d);...
         zeros(N,1)];
     
    x_z = reshape(x0,[n*N,1]);
    u_z = reshape(u0,[N,1]);
    
    z0 = [x_z;...
          u_z];
        
        
    for j=1:N-1
        r1 = (j-1)*n + 1;
        r2 = j*n; 

        Adiag{1,j} = A_d(r1:r2,:);
        Bdiag{1,j} = B_d(r1:r2,:);
        
        Idiag{1,j} = eye(n);
        
        F(r1:r2,1) = fx_d(r1:r2,1) - dt*A_c(r1:r2,:)*x0(:,j) - B_d(r1:r2,1)*u0(j);
    end
    
    Ishift = zeros(size(A_d));
    Ivals = blkdiag(Idiag{:});
    %Ivals = eye((N-1)*n);
    Imatrix = [Ishift,Ivals];
    
    M1 = [-blkdiag(Adiag{:}),Ishift] + Imatrix;
    
    M2 = -blkdiag(Bdiag{:});
    
    M = [M1,M2];

%end
            
% % Set up symbols for substitution
% x = sym('x',[n,1]);
% 
% % Loop through all temporal nodes
% for j=1:N-1
%     % Get proper indices for stacking A matrices at each node
%     r1 = (j-1)*n + 1;
%     r2 = j*n;
% 
%     % Verify divide by zero error
%     
%     % Replace divide by zero errors with eps for J
%     [Jnum,Jdenom] = numden(J_sym);
%     Jdiv0 = (double(subs(Jdenom,x,x0(:,j)))==0.0);
% 
%     if Jdiv0
%          fprintf('Divide by zero error in cost, replacing with eps!\n' )
%          J_sym = Jnum / eps;
%     end
%     
%     for i1=1:n
%         % Replace divide by zero errors with eps for dJ
%         [dJnum,dJdenom] = numden(dJ_sym(i1));
%         dJdiv0 = (double(subs(dJdenom,x,x0(:,j)))==0.0);
%         if dJdiv0
%             fprintf('Divide by zero error in cost deriv, replacing with eps!\n' )
%             dJ_sym(i1) = dJnum / eps
%         end
%     end
%     
%     % Determine continuous time numerical A matrix
%     J(r1:r2) = double(subs(J_sym,x,x0(:,j)));
% 
% 
%     % Determine continuous time numerical B matrix
%     dJ(r1:r2) = double(subs(dJ_sym,x,x0(:,j)));
% end
%}