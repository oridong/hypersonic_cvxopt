% Skye Mceowen
% ONR Update Simulations - Wang
% Jan15, 2021

clear all, close all, clc

fig_h = figure;

sigma_sweep = 0:20:180;
tf = 1000; % [s], final time

for ind=1:length(sigma_sweep)
    % Initialize control input values
        sigma = deg2rad(sigma_sweep(ind));  % [rad], bank angle
        T = 0;      % [N], thrust

    % Initialize state values (at burnout)
        h0 = 121e3; %98000;  % [m], initial altitude
        R = 6378e3; % [m], radius of the earth

        x0 = [R+h0;...          r [m]
             deg2rad(90);...    theta [rad]
             deg2rad(0);...     phi [rad]
             3700;...          V [m/s]
             deg2rad(16.4);...   gamma [rad]
             deg2rad(90);...    psi [rad]
             1200];           % m [kg]

    % Propogate dynamics
        [t,state_vec] = ode45(@(t,x) dynamics_wang(t,x,sigma,T), [0 tf], x0);

        t = t'; state_vec = state_vec';

    % Pull out variables
        r_vec       = state_vec(1,:)/1000;
        theta_vec   = state_vec(2,:);
        phi_vec     = state_vec(3,:);
        V_vec       = state_vec(4,:)/1000;
        gamma_vec   = state_vec(5,:);
        psi_vec     = state_vec(6,:);
        m_vec       = state_vec(7,:);

        sigma_vec = rad2deg(sigma)*ones(length(t),1);

    % Convert radius vector to cartesian value
        for i=1:length(r_vec)
            er(:,i) = e_r(theta_vec(i));
            r_xy(:,i) = r_vec(i)*er(:,i);
        end

    % Plot values
        figure(fig_h.Number)
        hold all
        circle(0,0,R/1000);
        hold all
        plot(r_xy(1,1),r_xy(2,1), "og")
        plot(r_xy(1,end),r_xy(2,end), "or")
        plot(r_xy(1,:),r_xy(2,:))
        title('Vehicle Approach to Earth')
        xlabel('x [km]')
        ylabel('y [km]')
        %xlim([-R 0]) 
        %ylim([-R inf])
        legend('Earths surface')

        figure(fig_h.Number+1);
        subplot(2,3,[1 2])
        hold all
        plot(t,r_vec-R/1000)
        title('Vehicle Altitude vs. Time')
        xlabel('Time [s]')
        ylabel('Altitude [km]')
        xlim([-inf inf])

        subplot(2,3,3)
        hold all
        plot(t,V_vec)
        title('Vehicle Velocity vs. Time')
        xlabel('Time [s]')
        ylabel('Velocity [km/s]')
        xlim([-inf inf])

        subplot(2,3,4)
        hold all
        plot(t,sigma_vec)
        ylim([-1,181])
        title('Control Input vs. Time')
        xlabel('Time [s]')
        ylabel('Sigma [deg]')
        xlim([-inf inf])


        subplot(2,3,5)
        hold all
        plot(t,rad2deg(theta_vec))
        title('Vehicle Theta vs. Time')
        xlabel('Time [s]')
        ylabel('Longitude [deg]')
        xlim([-inf inf])

        subplot(2,3,6)
        hold all
        plot(t,rad2deg(gamma_vec))
        title('Vehicle FPA vs. Time')
        xlabel('Time [s]')
        ylabel('Flight Path Angle [deg]')
        xlim([-inf inf])
    
end

figure(fig_h.Number)   
hold all
circle(0,0,R/1000)
legend('\sigma=0','\sigma=20','\sigma=40','\sigma=60','\sigma=80',...
        '\sigma=100','\sigma=120','\sigma=140','\sigma=160','\sigma=180','Earth')

figure(fig_h.Number+1)
subplot(2,3,4)
hold all
legend('\sigma=0','\sigma=20','\sigma=40','\sigma=60','\sigma=80',...
        '\sigma=100','\sigma=120','\sigma=140','\sigma=160','\sigma=180')

    
%% Functions
% Radial unit vector in XY plane (for plotting)
function e_r = e_r(theta)
    e_r = [ cos(theta) ;...
            sin(theta)];
end % end e_r

% Velocity unit vector in XY plane (for plotting)
function e_v = e_v(theta,fpa)
    phi_star = deg2rad(90) + fpa - theta;
    e_v = [ cos(phi_star) ;...
            sin(phi_star)];
end % end e_v

    

    
