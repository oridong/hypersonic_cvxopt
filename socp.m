% Skye Mceowen
% Qualifying Exam SOCP subproblem solver
% Hypersonic Reentry Vehicle Optimal Traj Gen
% Nov 6, 2020

function [O] = socp(I,ev,k_sc)    
    
    % Extract values from inputs struct
    dt = I.dt; 
    tf = I.tf;
    N = double(I.N);
    n = I.n;
    nN = double(N*(n+1));
    R = I.R;
    r_scale = ev.params.r_sf;
    
    u_min = I.u_min;
    u_max = I.u_max;
    
    f1max = I.f1max;
    f2max = I.f2max;
    f3max = I.f3max;
    
    x_i = I.x_i;
    x_f = I.x0(:,N);

    x0 = I.x0;  
    u0 = I.u0;
    
    % Pull out linearized matrices and cost, discretize
    [A_c, B_c, fx_c] = ev.linsys_c(ev.opt_in.A_hdl,ev.opt_in.B_hdl,x0,u0,N,n);    
    [A_d, B_d, fx_d] = ev.linsys_d(A_c,B_c,fx_c);

    [J_c, dJ_c] = ev.cost_c(ev.opt_in.J_hdl,ev.opt_in.dJ_hdl,x0);
    [J_d, dJ_d] = ev.cost_d(J_c,dJ_c);
    
    [f1,df1,f2,df2,f3,df3] = ev.path_constr(x0);
    
    fprintf('Finished with linearization!\n\n')
    
    % Pull out restacked cost and dynamic constraints!
    [c,z0,M,F] = ev.restack(x0,u0,dJ_d,A_d,B_d,fx_d,A_c,dt,N,n);
    
    fprintf('Finished with restack!\n\n')
    
%% BEGIN CVX
    cvx_tic
    cvx_begin
    % cvx_solver(I.cvx_solver)
    cvx_precision(I.cvx_precision)
    cvx_quiet(I.cvx_quiet)
    
    % Solution variables
    %variable z(nN,1)
    variable x(n,N)
    variable u(1,N)
    variable cost(1,N)

    
    % Minimization problem
    % TODO: Add full linearized cost
    % minimize(real(transpose(c)*z))
    % minimize(sum(cost))
    minimize(0)
    
    %for j=1:N
    %    r1 = (j-1)*n + 1;
    %    r2 = j*n;
    %    cost(j) == real(transpose(c(r1:r2)))*x(:,j);
    %end
    
    % Initial conditions
    %z(1,1) <= x_i(1);
    %1 <= z(1,1);
    x(:,1) == x_i;
    
    
    % Final conditions
    % TODO: Select sensible final conditions
    %x(:,N) == x_f;
    
    % Dynamics
    % M*z == F;
    %%{
    for j=N
        %r1 = (j-1)*n + 1;
        %r2 = j*n;
        %j_u = n*N + j;
        
        if j<N
        x(:,j+1) - A_d(r1:r2,:)*x(:,j) - B_d(r1:r2,1)*u(:,j) ==  ...
                       fx_c(r1:r2) - A_d(r1:r2,:)*x0(:,j) - B_d(r1:r2,1)*u0(:,j);
        end
        % Verifying cost is formed via states
        %z(r1:r2,1) == x(:,j);          
        %z(j_u,1) == u(:,j);
        
        % State constraints
        %norm(x(1,j)) <= ev.ic.r_i;
        %1 <= x(1,j);
        %norm(x(3,j)) <= 1.2;
    end
    %}
   
    
    % State constraints
    %for jz=1:nN
        %norm(z(jz)) <= z0(1);
    %end
    
    % Control constraints
    for j=1:N
        %j_u = n*N + j;
        %norm(z(j_u,1),2) <= 1 ;%u_max;
    %    u_min <= u(1,j) <= u_max;
        norm(u(1,j)) <= 1
    end
    
    % Path constraints
    %for j=1:N
    %    f1 + df1(j,1)*(x(1,j)-x0(1,j)) + df1(j,2)*(x(3,j)-x0(3,j)) <= f1max;
    %    f2 + df2(j,1)*(x(1,j)-x0(1,j)) + df2(j,2)*(x(3,j)-x0(3,j)) <= f2max;
    %    f3 + df3(j,1)*(x(1,j)-x0(1,j)) + df3(j,2)*(x(3,j)-x0(3,j)) <= f3max;
    %end
    
    % Trust region constraints
    for j=1:N
        I.eps_conv(3) = abs(1/x0(3,j));
        for ind = 1:n
           norm(x(n,j) - x0(n,j)) <= I.delta_tr(n); 
        end
    end
    
    
    cvx_end
    %
    % END CVX
    %
    
    % Pull out x from z
    %x = z(1:n*N);
    %x = reshape(x,[n,N]);
    %u = z(n*N+1:end)';
    
    
    % Convergence check
    for i=1:n
        for j=1:N 
            conv_check(i,j) = norm(x(i,j) - x0(i,j)); 
        end 
        [conv_max(i,1) idx(i,1)] = max(conv_check(i,:));
    end
    I.eps_conv(3) = abs(500/x0(3,idx(3)));
    converged_bool = (conv_max <=  I.eps_conv);
    
    
    % Output Assignment
    O = struct;
    O.t = linspace(0,tf,N);
    O.x = full(x);
    O.u = full(u);
    %O.z = full(z);
    O.cost = full(cost); % transpose(c)*z;
    O.converged = min(converged_bool); %O.converged = min(min(O.x - x0 < I.eps_conv));
    O.status = cvx_status;
    O.timing = cvx_toc;
    
    % Print statements
    fprintf('iter = %02d, ',k_sc)
    fprintf('conv = %01d, ',O.converged)
    fprintf('tsol = %04d, ',ceil(1000*O.timing(5)))
    fprintf('status = %s',cvx_status)
    fprintf('\n')
end