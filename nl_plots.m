% Skye Mceowen
% Qualifying Exam Plots
% Nov4, 2020

clear all, close all, clc
% Select directory (later use if statement to pick quals vs wang)
addpath('matfiles/')
% addpath('qual_fcns/')


% Created fixed bank angle plots

sigma_vec = 0:20:180;
h = figure;
for ind=1:length(sigma_vec)
    %
    % Select bank angle sigma in degrees!
    %
    sigma0 = deg2rad(sigma_vec(ind));      % [rad], initial control
    
    %
    % Define and initialize values
    %
    % Create vehicle
    ev = vehicle(sigma0);

    % Problem parameters
    ev.opt_in.dt = 0.8*ev.params.t_sf; 
    ev.opt_in.tf = 2500*ev.params.t_sf;             % [s], final time
    ev.opt_in.N = round(ev.opt_in.tf/ev.opt_in.dt +1);     % [-], number of discrete points


    %
    % Propogate initial conditions with full NL dynamics
    %   to produce initial trajectory for solver
    % 
    % Initialize trajectory
    [t,r0,theta0,v0,fpa0] = ev.gen_traj();

    % Create initial x vector
    t_vec{ind} = t;
    x0{ind} = [r0;theta0;v0;fpa0];
    u0{ind} = ev.ic.u_i*ones(1,ev.opt_in.N);

    ev.opt_in.x0 = x0{ind};
    ev.opt_in.u0 = u0{ind};

    ev.plot_traj(t,x0{ind},u0{ind},h) 
end

subplot(2,3,4)
hold all
legend('\sigma=0','\sigma=20','\sigma=40','\sigma=60','\sigma=80',...
        '\sigma=100','\sigma=120','\sigma=140','\sigma=160','\sigma=180')
    
savefig('matfiles/BankAngle.fig')
save('matfiles/bankangle_figs.mat')