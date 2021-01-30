% Skye Mceowen
% Qualifying Exam Vehicle Class Defition
% Hypersonic Reentry
% Nov6, 2020

classdef vehicle
    %
    % Member variables
    %
    properties
        % Values
        %Value {mustBeNumeric}
        params = struct; 
        ic = struct;
        opt_in = struct;
        opt_out = struct;
        sysdyn = vehicle_dynamics;
        
        % Functions
        fn = struct;
    end
    

    %
    % Member functions
    %
    methods
%% CLASS INITIALIZATION
        %
        % Constructor
        %
        function this = vehicle(sigma_i)
            if nargin==0
               sigma_i = deg2rad(180);
            end
            %
            % Constants
            %
            %Physical constants
            this.params.R =  6371; % [km], radius of earth
            this.params.M = 5.972e24; % [kg], mass of the earth
            this.params.g0 = 9.81/1000; % [km/s^2], gravity at earth surface;
            this.params.G = 6.67408e-11/(1000)^3; % [km^3/(kg s^2)], gravitational constant
            this.params.rho0 = 1.3*(1000)^3; %kg/km^3 density
            this.params.H = 7; % k*T/(mbar*g), approx 7000m up to 100km
            this.params.mu = this.params.G*this.params.M; %km gravitational standard param

            % Vehicle parameters
            this.params.A = 291.22 / (1000)^2; % [km^2], reference area
            this.params.m = 104.305; % [kg], mass 
            this.params.r_sf = 1/this.params.R;
            this.params.v_sf = 1/sqrt(this.params.R*this.params.g0);
            this.params.t_sf = 1/sqrt(this.params.R/this.params.g0);

            % Problem parameters
            this.ic.h_i = 121.9; % [km], initial altitude
            r_i = this.ic.h_i + this.params.R; % [km], initial radius
            this.ic.r_i = r_i * this.params.r_sf;
            this.ic.theta_i = deg2rad(90); % [rad], angle of initial r_i in CCW from X axis

            v_i = 7.627; % [km/s], initial velocity magnitude
            this.ic.v_i = v_i * this.params.v_sf;
            this.ic.fpa_i = deg2rad(-5);%deg2rad(-0.5); % [rad], FPA of initial v_i 

            this.ic.x_i = [ this.ic.r_i;...     x1
                            this.ic.theta_i;... x2
                            this.ic.v_i;...     x3
                            this.ic.fpa_i];     %x4
            
            this.ic.u_i = cos(sigma_i);
            %%
            %
            % Function handles
            %
            %TODO: Fix piecewise things later
            
            
            % RHO
            this.fn.h = @(r) r/this.params.r_sf - this.params.R;
            this.fn.rho_hdl = @(r) this.params.rho0*exp(-this.fn.h(r)/this.params.H);
                            %{
                            = @(r) piecewise(r > 1, ... if
                                           this.params.rho0*exp(-this.fn.h(r)/this.params.H),...
                                        r <= 1, ... elseif
                                           this.params.rho0);
                            %}
                                        
                                        
            % ALPHA                    
            this.fn.alpha_hdl = @(v) 40 - 0.20705*(1000*v-4570)^2/(340^2); 
                            %{
                            = @(v) piecewise(v>4.570, ... if
                                                40,... then
                                            v<=4.570, ... elseif
                                                40-0.20705*(1000*v-4570)^2/(340^2)); % then  
                            %}
            
            
            % LIFT COEFFICIENT    
            this.fn.Cl_hdl = @(v) -0.041065 ...
                                    + 0.016292*this.alpha(v) ...
                                    + 0.0002602*(this.alpha(v)^2);
            % DRAG COEFFICIENT
            this.fn.Cd_hdl = @(v) 0.080505 ...
                                    - 0.03026*this.fn.Cl_hdl(v) ...
                                    + 0.86495*(this.fn.Cl_hdl(v))^2;

            % LIFT
            m = this.params.m;
            R = this.params.R;
            A = this.params.A;
            g0 = this.params.g0;
            H = this.params.H;
            this.fn.L_hdl = @(r,v) (A*R/(2*m)) * this.rho(r) ...
                                                * v^2 * this.fn.Cl_hdl(v);

            % DRAG
            this.fn.D_hdl = @(r,v) (A*R/(2*m)) * this.rho(r) ...
                                                * v^2 * this.fn.Cd_hdl(v);


            % NONLINEAR DYNAMICS
            % NOTE : x1=r, x2=theta, x3=v, x4=gamma
            this.fn.fx_hdl = @(x,u) ... 
                            [   x(3) * sin(x(4));... x1Dot
                                x(3) * cos(x(4))/(x(1));... x2Dot
                                -this.fn.D_hdl(x(1),x(3)) - sin(x(4))/(x(1)^2) ;...x3Dot
                                ( x(3)/x(1) - 1/((x(1)^2)*x(3)) ) * cos(x(4)) ...  x4Dot
                                    + u*this.fn.L_hdl(x(1),x(3))/x(3)];

            % COST FUNCTION
            % Min heat rate
            % NOTE : x1=r, x2=theta, x3=v, x4=gamma
            kq = 9.4369e-5;
            sqrtR0g0 = 1/this.params.v_sf;
            k0 = kq*sqrtR0g0^(3.15);    
            k0_inv = kq/sqrtR0g0^(3.15);  
            this.fn.cost_hdl = @(x) k0*sqrt(this.fn.rho_hdl(x(1))) * x(3)^3.15;
            
            
            % PATH CONSTRAINTS
            % NOTE: this.params.v_sf = 1/sqrt(this.params.R*this.params.g0);
            % Heat rate
            this.fn.f1_hdl = @(r,v) kq*sqrt(this.rho(r)) * v^(3.15);
            this.fn.df1_hdl = @(r,v)  [ -0.5*R*k0_inv*v^(3.15) * sqrt(this.rho(r))/H,...
                                        3.15*k0*v^(2.15) * sqrt(this.rho(r))];

            % Dynamic pressure
            this.fn.f2_hdl = @(r,v) 0.5 * this.rho(r) * v^2;
            this.fn.df2_hdl = @(r,v)  [ -0.5*R*this.rho(r)*(v*sqrtR0g0)^2 / H, ...
                                        sqrtR0g0^2 * this.rho(r) * v];

            % Normal load
            this.fn.f3_hdl = @(r,v) sqrt(this.fn.D_hdl(r,v)^2 + this.fn.D_hdl(r,v)^2);
            this.fn.df3_hdl = @(r,v) [ -0.5 * R^2 * A * this.fn.f3_hdl(r,v) * this.rho(r) * v^2 /(m*H),...
                                        R * A * this.fn.f3_hdl(r,v) * this.rho(r) * v / m]; 


            %%
            %
            % Optimization inputs
            %
            this.opt_in.n = 4;
            this.opt_in.N = 100;
            this.opt_in.tf = 1600 * this.params.t_sf;
            this.opt_in.dt = this.opt_in.tf/(this.opt_in.N - 1);


            [Ahdl, Bhdl] = linsys_sym(this);
            this.opt_in.A_hdl = Ahdl;
            this.opt_in.B_hdl = Bhdl;
            
            % TODO: Add cost fcn values
            [Jhdl, dJhdl] = cost_sym(this);
            this.opt_in.J_hdl = Jhdl;
            this.opt_in.dJ_hdl = dJhdl;
            
            this.opt_in.x0 = this.ic.x_i;
            this.opt_in.u0 = this.ic.u_i;
            
            % CVX params
            % TODO: Install ECOS
            %this.opt_in.cvx_solver = 'ecos';             % solver
            this.opt_in.cvx_precision = 'low';  % precision
            this.opt_in.cvx_quiet = true;       % cvx print option
            
            % Convergence params
            this.opt_in.k_max = 20;                     % max number of successive iterations
            this.opt_in.lcvx_tol = 1e-2;                % losslessness tolerance [N]
            this.opt_in.eps_conv = [ 0.1/R;...  % convergence tolerance
                                     deg2rad(0.05);...
                                     0.001/(this.ic.v_i/this.params.v_sf);...
                                     deg2rad(0.05)];
            
            % Trust region params
            this.opt_in.delta_tr = [ 10/R;...  % trust region convergence tolerance
                                     deg2rad(20);...
                                     5/(this.ic.v_i/this.params.v_sf);...
                                     deg2rad(20)];
            
            
            % Physical params
            this.opt_in.m = this.params.m;
            this.opt_in.R = this.params.R;
            this.opt_in.g0 = this.params.g0;
            
            % Control input bounds
            this.opt_in.u_min = -1;
            this.opt_in.u_max =  1;
            
            % Boundary conditions
            % TODO: Fix this terminal BC
            this.opt_in.x_i = this.ic.x_i; % Initial BC
            this.opt_in.x_f = 0; % terminal BC
            
            % Constraint limits
            this.opt_in.f1max = 1500*(1000)^2;  % [kW/km^2], max heat rate
            this.opt_in.f2max = 18,000*(1000)^2;% [N/km^2], max q 
            this.opt_in.f3max = 2.5*g0; % [km/s^], max load

        end % end constructor
      
      
      
        % Resetting initial condition
        function reset_ic(this,x_i)
            if nargin==1
                x_i = [  121.9 + this.params.R;...
                        deg2rad(90);...
                        7.627;...
                        deg2rad(-2)];
            end

            this.ic.r_i = x_i(1)/this.params.r_sf;
            this.ic.theta_i = x_i(2);

            this.ic.v_i = x_i(3)/this.params.v_sf;           
            this.ic.aoa_i = x_i(4);
        end % end reset_ic
        
        
%% DYNAMICS and COST FUNCTIONS     

        function rho = rho(this,r)
            % Check if numeric
            if isnumeric(r)                                     
                if norm(r)>1
                    rho = this.fn.rho_hdl(r);
                else
                    rho = this.params.rho0;
                end
           
            % Otherwise, if symbolic:
            else
               rho = this.fn.rho_hdl(r); 
            end
        end % end fx

        function alpha = alpha(this,v)
           % Check if numeric
            if isnumeric(v)                                     
                if norm(v)>4.570
                    alpha = 40;
                else
                    alpha = this.fn.alpha_hdl(v);
                end
           
            % Otherwise, if symbolic:
            else
                alpha = this.fn.alpha_hdl(v); 
            end
        end
        
        function fx = fx(this,x0,u0)
            fx = this.fn.fx_hdl(x0,u0);
            
            for j=1:length(x0)
                if isinf(fx(j));
                   fx(j) = 1/eps; 
                end
            end
        end % end fx
        
        
        % Determine symbolic linearized A matrix
        function [A_hdl, B_hdl] = linsys_sym(this)
            % Create symbolic values
            n = this.opt_in.n;
            syms u
            x = sym('x',[n 1]);
           
            % Create symbolic derivatives vector
            derivs = this.fn.fx_hdl(x,u);
           
            % Compute Jacobian matrices
            A_sym = jacobian(derivs,x);
            B_sym = jacobian(derivs,u);
            
            A_hdl = matlabFunction(A_sym,'Vars', {[x; u]});
            B_hdl = matlabFunction(B_sym,'Vars', {[x; u]});
        end % end linsys_sym
        
        
        
        function [A_c, B_c, fx_c] = linsys_c(this,A_hdl, B_hdl,...
                                                x0,u0,N,n)
                                    
            % Set number of temporal nodes N if not given
            if nargin==5
               [n, N] = size(x0);
               opt_in.N = N; % temporal nodes
               opt_in.n = n; % length of state vector
            end 

            
            % Loop through all temporal nodes
            for j=1:N
                % Get proper indices for stacking A matrices at each node
                r1 = (j-1)*n + 1;
                r2 = j*n;
                
                % Determine continuous time numerical fx
                fx_c(r1:r2,:) = this.fx(x0(:,j),u0(j));
                
                
                % Determine continuous time numerical A matrix
                Amat = A_hdl([x0(:,j); u0(j)]);
                
                
                % Determine continuous time numerical B matrix
                Bmat = B_hdl([x0(:,j); u0(j)]);
                
                
                % Remove Inf and NaN values
                for row=1:n
                    if isinf(Bmat(row))
                        Bmat(row) = 1 / eps;
                    elseif isnan(Bmat(row))
                        Bmat(row) = 0;
                    end 
                    for col=1:n 
                        if isinf(Amat(row,col))
                            Amat(row,col) = 1 / eps;
                        elseif isnan(Amat(row,col))
                            Amat(row,col) = 0;
                        end 
                    end
                end
                
                
                % Set matrices
                A_c(r1:r2,:) = Amat;
                B_c(r1:r2,1) = Bmat;
            end % end for j=1:N
            
        end
        
        % Numerical A,B values for discrete time dynamics
        function [A_d, B_d, fx_d] = linsys_d(this,A_c,B_c,fx_c,dt)
            if nargin==4
               dt = this.opt_in.dt; 
            end
            
            % Pull out size of state space
            [m, n] = size(A_c);

            % Pull out number of temporal nodes
            %N = m/n + 1;
            N = m/n;
            
            % Create discrete fx_d
            fx_d = dt*fx_c;
                                    
            for j=1:N
                % Get proper indices for stacking A matrices at each node
                r1 = (j-1)*n + 1;
                r2 = j*n;               
                
                
                % Create state space model for continuous linear system
                cont_sys = ss(A_c(r1:r2,:),B_c(r1:r2),[],[]);
                
                                
                % Convert to discrete-time model
                %dis_sys = c2d(cont_sys,dt);
                %A_d(r1:r2,:) = double(dis_sys.A);
                %B_d(r1:r2,1) =  double(dis_sys.B);
                
                
                % TESTING
                % TODO: When breaking here and exponentialting Ac^n
                %   discovered A(3,3) and A(4,3) blow up to infty
                %   which indicates velocity and gamma coupling
                %   is unstable ... look at dynamics!
                %   (could also be due to u term...)
                A_d(r1:r2,:) = eye(n) + dt* double(cont_sys.A );
                B_d(r1:r2,1) =  dt* double(cont_sys.B );
            end
            
        end % end linsys_d
        
        
        % Determine symbolic linearized cost
        function [cost_hdl, dcost_hdl] = cost_sym(this)
            % Create symbolic values
            n = this.opt_in.n;
            x = sym('x',[n 1]);
           
            % Create symbolic derivatives vector
            cost_sym = this.fn.cost_hdl(x);
            cost_hdl = this.fn.cost_hdl;
           
            % Compute Jacobian matrices
            dcost_sym = jacobian(cost_sym,x);
            dcost_hdl = matlabFunction(dcost_sym, 'Vars', {x});
        end % end cost_sym
        
        % Numerical continuous cost
        function [J_c, dJ_c] = cost_c(this,J_hdl,dJ_hdl,x0,N,n)
            % Set number of temporal nodes N if not given
            if nargin==4
               [n, N] = size(x0);
               opt_in.N = N; % temporal nodes
               opt_in.n = n; % length of state vector
            end 
            
            
            for j=1:N
                % Get proper indices for stacking A matrices at each node
                r1 = (j-1)*n + 1;
                r2 = j*n;


                % Determine continuous time numerical A matrix
                J_c(j,1) = J_hdl(x0(:,j));
                
                
                % Remove NaNs and inf
                if isnan(J_c(j))
                    J_c(j) = 0;
                elseif isinf(J_c(j))
                    J_c(j) = 1 / eps;
                end


                % Determine continuous time numerical B matrix
                dJ_c(1,r1:r2) = dJ_hdl(x0(:,j));
            end
        end
        
        % Numerical discrete cost
        function [J_d, dJ_d] = cost_d(this,J_c,dJ_c,dt)
            if nargin==3
               dt = this.opt_in.dt; 
            end
            
            J_d = dt*J_c;
            dJ_d = dt*dJ_c;
        end
        
        
        function [f1,df1,f2,df2,f3,df3] = path_constr(this,x0)
            [n,N] = size(x0);
            
            for j = 1:N
                % Pull out r and v
                r = x0(1,j);
                v = x0(3,j);
                
                %
                % Heat rate
                %
                f1(j,1) = this.fn.f1_hdl(r,v);
                df1(j,:) = this.fn.df1_hdl(r,v);
                
                
                if isinf(f1(j,1))
                    f1(j,1) = sign(f1(j,1)) / eps;
                end
                
                if isinf(df1(j,1))
                    df1(j,1) = sign(df1(j,1)) / eps;
                end
                
                if isinf(df1(j,2))
                    df1(j,2) = sign(df1(j,2)) / eps;
                end

                %
                % Dynamic pressure
                %
                f2(j,1) = this.fn.f2_hdl(r,v);
                df2(j,:) = this.fn.df2_hdl(r,v);
                
                
                if isinf(f2(j,1))
                    f2(j,1) = sign(f2(j,1)) / eps;
                end
                
                if isinf(df2(j,1))
                    df2(j,1) = sign(df2(j,1)) / eps;
                end
                
                if isinf(df2(j,2))
                    df2(j,2) = sign(df2(j,2)) / eps;
                end
                

                %
                % Normal load
                %
                f3(j,1) = this.fn.f3_hdl(r,v);
                df3(j,:) = this.fn.df3_hdl(r,v);
                
                
                if isinf(f3(j,1))
                    f3(j,1) = sign(f3(j,1)) / eps;
                end
                
                if isinf(df3(j,1))
                    df3(j,1) = sign(df3(j,1)) / eps;
                end
                
                if isinf(df3(j,2))
                    df3(j,2) = sign(df3(j,2)) / eps;
                end
                
                f1 = real(f1); 
                df1 = real(df1);
                f2 = real(f2); 
                df2 = real(df2);
                f3 = real(f3);
                df3 = real(df3);
            end    
        end
        
        
        % Reorganize cost and constraints into 
        function [c,z0,M,F] = restack(this,x0,u0,dJ_d,A_d,B_d,fx_d,A_c,dt,N,n)
            if nargin==7
               [n, N] = size(x0);
               dt = opt_in.dt;
            end 
            
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
        end
        
%% SIMULATION
        function xDot = derivs(this,t,x,u)
            xDot = this.fx(x,u);
        end

        % Generating trajectories through NL simulation
        % TODO: Allow input arguments for arbitrary ICs
        function [tj,r0,theta0,v0,fpa0] = gen_traj(this)
            % Initialize inputs
            u_i = this.ic.u_i;
            x_i = this.ic.x_i;
            dt  = this.opt_in.dt;
            tf  = this.opt_in.tf;
            derivs_hdl = @(t,x) this.derivs(t,x,u_i);
            
            % Propogate state dynamics
            [tj,state_vec] = ode45(@(t,x) derivs_hdl(t,x), 0:dt:tf, x_i);

            tj = tj'; state_vec = state_vec';
            
            this.opt_in.times = tj;
            this.opt_in.x0 = state_vec;
            
            r0      = state_vec(1,:);
            theta0  = state_vec(2,:);
            v0      = state_vec(3,:);
            fpa0    = state_vec(4,:);
        end % end gen_traj
        
        % Radial unit vector in XY plane (for plotting)
        function e_r = e_r(this,theta)
            e_r = [ cos(theta) ;...
                    sin(theta)];
        end % end e_r
        
        % Velocity unit vector in XY plane (for plotting)
        function e_v = e_v(this,theta,fpa)
            phi_star = deg2rad(90) + fpa - theta;
            e_v = [ cos(phi_star) ;...
                    sin(phi_star)];
        end % end e_v
        
        % Plotting
        % TODO: Allow input arguments for arbitrary trajectories
        function plot_traj(this,t,x0,u0,fig_h)
            if nargin == 1
                % Initialize inputs
                t     = this.opt_in.times;
                x0    = this.opt_in.x0;
                u0    = this.opt_in.u0;
                fig_h = figure;
            elseif nargin==3
                u0    = this.opt_in.u0;
                fig_h = figure;
            elseif nargin==4
                fig_h = figure;
            end
            
        
                r0      = x0(1,:)/this.params.r_sf;
                theta0  = x0(2,:);
                v0      = x0(3,:)/this.params.v_sf;
                fpa0    = x0(4,:);
                t = t/this.params.t_sf;
                
            for i=1:length(theta0)
                er(:,i) = this.e_r(theta0(i));
                r_vec(:,i) = r0(i)*er(:,i);
                
                ev(:,i) = this.e_v(theta0(i),fpa0(i));
                v_vec(:,i) = v0(i)*ev(:,i);
            end
            
            figure(fig_h.Number);
            % Plot result 
            subplot(2,3,1)
            hold all
            circle(0,0,this.params.R);
            hold all
            plot(r_vec(1,:),r_vec(2,:))
            title('Vehicle Approach to Earth')
            xlabel('x [km]')
            ylabel('y [km]')
            xlim([-this.params.R 0]) 
            ylim([-this.params.R inf])
            legend('Earths surface')

            subplot(2,3,2)
            hold all
            plot(t,r0-this.params.R)
            title('Vehicle Altitude vs. Time')
            xlabel('Time [s]')
            ylabel('Altitude [km]')
            xlim([-inf inf])

            subplot(2,3,3)
            hold all
            plot(t,v0)
            title('Vehicle Velocity vs. Time')
            xlabel('Time [s]')
            ylabel('Velocity [km/s]')
            xlim([-inf inf])

            subplot(2,3,4)
            hold all
            plot(t,acosd(u0))
            ylim([-1,181])
            title('Control Input vs. Time')
            xlabel('Time [s]')
            ylabel('Bank Angle [deg]')
            xlim([-inf inf])
            
            
            subplot(2,3,5)
            hold all
            plot(t,rad2deg(theta0))
            title('Vehicle Theta vs. Time')
            xlabel('Time [s]')
            ylabel('Longitude [deg]')
            xlim([-inf inf])

            subplot(2,3,6)
            hold all
            plot(t,rad2deg(fpa0))
            title('Vehicle FPA vs. Time')
            xlabel('Time [s]')
            ylabel('Flight Path Angle [deg]')
            xlim([-inf inf])
            
        end % end plot_traj
        
        function plot_bank_sweep(this,fig_h)
            % Created swept fixed bank angle plots
            sigma_vec = 0:20:180;
            for ind=1:length(sigma_vec)
                %
                % Select bank angle sigma in degrees!
                %
                sigma0 = deg2rad(sigma_vec(ind));      % [rad], initial control

                %
                % Define and initialize values
                %
                % Create vehicle
                this.ic.u_i = cos(sigma0);

                % Problem parameters
                this.opt_in.dt = 0.8*this.params.t_sf; 
                this.opt_in.tf = 2500*this.params.t_sf;             % [s], final time
                this.opt_in.N = round(this.opt_in.tf/this.opt_in.dt +1);     % [-], number of discrete points


                %
                % Propogate initial conditions with full NL dynamics
                %   to produce initial trajectory for solver
                % 
                % Initialize trajectory
                [t,r0,theta0,v0,fpa0] = this.gen_traj();

                % Create initial x vector
                t_vec{ind} = t;
                x0_vec{ind} = [r0;theta0;v0;fpa0];
                u0_vec{ind} = this.ic.u_i*ones(1,this.opt_in.N);

                this.opt_in.x0 = x0_vec{ind};
                this.opt_in.u0 = u0_vec{ind};

                this.plot_traj2(t,x0_vec{ind},u0_vec{ind},fig_h) 
            end
            
            figure(fig_h.Number)   
            hold all
            circle(0,0,this.params.R)
            legend('\sigma=0','\sigma=20','\sigma=40','\sigma=60','\sigma=80',...
                    '\sigma=100','\sigma=120','\sigma=140','\sigma=160','\sigma=180','Earth')

            figure(fig_h.Number+1)
            subplot(2,3,4)
            hold all
            legend('\sigma=0','\sigma=20','\sigma=40','\sigma=60','\sigma=80',...
                    '\sigma=100','\sigma=120','\sigma=140','\sigma=160','\sigma=180')

        end %end plot_bank_sweep

                % Plotting 2
        function plot_traj2(this,t,x0,u0,fig_h)
            if nargin == 1
                % Initialize inputs
                t     = this.opt_in.times;
                x0    = this.opt_in.x0;
                u0    = this.opt_in.u0;
            elseif nargin==3
                u0    = this.opt_in.u0;
            elseif nargin==4
                fig_h = figure;
            end


                r0      = x0(1,:)/this.params.r_sf;
                theta0  = x0(2,:);
                v0      = x0(3,:)/this.params.v_sf;
                fpa0    = x0(4,:);
                t = t/this.params.t_sf;

            for i=1:length(theta0)
                er(:,i) = this.e_r(theta0(i));
                r_vec(:,i) = r0(i)*er(:,i);

                ev(:,i) = this.e_v(theta0(i),fpa0(i));
                v_vec(:,i) = v0(i)*ev(:,i);
            end

            figure(fig_h.Number);
            % Plot result 
            hold all
            plot(r_vec(1,:),r_vec(2,:))
            title('Vehicle Approach to Earth')
            xlabel('x [km]')
            ylabel('y [km]')
            xlim([-this.params.R 0]) 
            ylim([-this.params.R inf])
            %xlim([-inf inf]) %xlim([min(r_vec(1,:)),max(r_vec(1,:))])
            %ylim([0.98*this.params.R,1.001*max(r_vec(2,:))]) %ylim([min(r_vec(2,:)),max(r_vec(2,:))])
            %legend('Earths surface')

            figure(fig_h.Number+1);
            subplot(2,3,[1 2])
            hold all
            plot(t,r0-this.params.R)
            title('Vehicle Altitude vs. Time')
            xlabel('Time [s]')
            ylabel('Altitude [km]')
            xlim([-inf inf])

            subplot(2,3,3)
            hold all
            plot(t,v0)
            title('Vehicle Velocity vs. Time')
            xlabel('Time [s]')
            ylabel('Velocity [km/s]')
            xlim([-inf inf])

            subplot(2,3,4)
            hold all
            plot(t,acosd(u0))
            ylim([-1,181])
            title('Control Input vs. Time')
            xlabel('Time [s]')
            ylabel('Bank Angle [deg]')
            xlim([-inf inf])


            subplot(2,3,5)
            hold all
            plot(t,rad2deg(theta0))
            title('Vehicle Theta vs. Time')
            xlabel('Time [s]')
            ylabel('Longitude [deg]')
            xlim([-inf inf])

            subplot(2,3,6)
            hold all
            plot(t,rad2deg(fpa0))
            title('Vehicle FPA vs. Time')
            xlabel('Time [s]')
            ylabel('Flight Path Angle [deg]')
            xlim([-inf inf])

        end % end plot_traj2
        
        
        % Inputs - TODO: Install ECOS
        function I = inputs(this,x0,u0)
            %if nargin==1
            %    x0 = this.opt_in.x0(1,:);
            %    u0 = this.opt_in.u0;
            %end
            
            % Boundary conditions
            % TODO: Fix this terminal BC
            this.opt_in.N = double(this.opt_in.N);
            this.opt_in.x_i = x0(:,1); % Initial BC
            
            % Initialization Trajectory
            this.opt_in.x0 = x0;
            this.opt_in.u0 = u0*ones(1,int64(this.opt_in.N));
            
            % Control constraints
            this.opt_in.u_max = cosd(0);
            this.opt_in.u_min = cosd(80);
            
            
            I = this.opt_in;
            
        end %end inputs
        
    end % end methods
end % end classdef

%
% Function format:
%
%{

% DEFINE FUNCTIONS IN CONSTRUCTOR AS: 
this.fn.FCN_hdl = @(INPUT) piecewise(INPUT >= , ),...
                        INPUT < , );  

% Check if numeric
if isnumeric(INPUT)
    syms INPUTsym
    FCN_sym = this.fn.FCN_hdl(INPUTsym);
    OUTPUT = double(subs(FCN_sym,INPUTsym,INPUT));

% Otherwise, if symbolic:
else
    OUTPUT = this.fn.FCN_hdl(INPUT); 
end

%}