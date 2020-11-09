% Skye Mceowen
% Qualifying Exam Simulations
% Nov4, 2020

clear all, close all, clc
% Select directory (later use if statement to pick quals vs wang)
addpath('matfiles/')
% addpath('qual_fcns/')

%
% Select bank angle sigma in degrees!
%
sigma0 = deg2rad(180);      % [rad], initial control

% Create vehicle
ev = vehicle(sigma0);

%
% Define and initialize values
%
dt = 0.08*ev.params.t_sf;

N = 1000;
tf = dt*(N-1);

%tf = 100*ev.params.t_sf;             % [s], final time
%N = round(tf/dt +1);     % [-], number of discrete points

% LOOKS GOOD
%dt = 0.8*ev.params.t_sf; 
%tf = 2500*ev.params.t_sf;             % [s], final time
%N = round(tf/dt +1);     % [-], number of discrete points

% Problem parameters
ev.opt_in.dt = dt; 
ev.opt_in.N = N; 
ev.opt_in.tf = tf;   

%
% Propogate initial conditions with full NL dynamics
%   to produce initial trajectory for solver
% 

% Initialize trajectory
[t,r0,theta0,v0,fpa0] = ev.gen_traj();
% save('x0_traj3_N100.mat','t','r0','theta0','v0','fpa0')
%load('x0_traj1600s.mat')

% Create initial x vector
x0 = [r0;theta0;v0;fpa0];
u0 = ev.ic.u_i;

ev.opt_in.x0 = x0;
ev.opt_in.u0 = u0;

ev.plot_traj(t,x0,ones(1,N)*u0)

%%{
%
% Solve SOCP problems until convergence
%
% Initialize inputs
I = ev.inputs(x0,u0);
I = {I};
O = cell(0,1);
solve_time = 0;


% Sequential SOCP Loop
for k=1:I{1}.k_max
    O{k,1} = socp(I{k,1},ev,k);
    solve_time = solve_time+O{k}.timing(5);
  
    if (O{k}.converged) 
    % Check convergence.
    %   If so, z^* = z_k
    fprintf('converged, total solve time = %04d [ms]\n', ...
        ceil(1000*solve_time))
        O{k}.x = full(O{k}.x);
        fprintf('\n')
        break
    elseif (length(O{k}.status)) == length('Infeasible') % == 'Infeasible')
        if (O{k}.status == 'Infeasible')
            fprintf('Sorry, problem infeasible\n!')
            break 
        end
    elseif (k < I{1}.k_max)
        if (length(O{k}.status) == length('Failed')) % == 'Failed')
            if (O{k}.status == 'Failed')
                fprintf('Sorry, optimizer failed!\n')
         %       break 
            end
        end
        % If not, k = k+1 and repeat optimization.
        I{k+1,1} = I{k};
        I{k+1,1}.x0 = O{k,1}.x;
        I{k+1,1}.u0 = O{k,1}.u;
        I{k+1,1}.k = k;
  end
end



lng = length(O);
ind = lng - 1;
ev.plot_traj(O{ind}.t,O{ind}.x,O{ind}.u)

save('matfiles/SUCCESS_9.mat')
%}

%
% TODO !
%
% TODO: MAKE ALL INITIAL STATES UNITLESS?
% TODO: Scale EVERYTHING properly  
% TODO: Verify scaling in cost function...
% TODO: After scaling, relinearize cost
% TODO: Stack cost into c^Tz ... 
% TODO: Go through and add other constraints
% TODO: Simulate for various initial sigma0 values
% TODO: Possibly 
% TODO: Use all physical values from paper


