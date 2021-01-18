% Skye Mceowen
% ONR Update Simulations
% Jan15, 2021

clear all, close all, clc

fig_h = figure;

alpha_sweep = 0:20:180;
tf = 1000; % [s], final time

for ind=1:length(alpha_sweep)
    % Initialize control input values
        alpha = deg2rad(alpha_sweep(ind));  % [rad], angle of attack
        sigma = deg2rad(0);  % [rad], bank angle
        T = 0;      % [N], thrust

    % Initialize state values (at burnout)
        h0 = 121e3; %98000;  % [m], initial altitude
        R = 6378e3; % [m], radius of the earth

        x0 = [R+h0;...          r [m]
             deg2rad(90);...    theta [rad]
             deg2rad(0);...     phi [rad]
             3700;...           V [m/s]
             deg2rad(16.4);...  gamma [rad]
             deg2rad(90);...    psi [rad]
             1200;...           m [kg]
             0];                % Q_init

    % Propogate dynamics
        [t,state_vec] = ode45(@(t,x) dynamics(t,x,sigma,alpha,T), [0 tf], x0);

        t = t'; state_vec = state_vec';

    % Pull out variables
        r_vec       = state_vec(1,:)/1000;
        theta_vec   = state_vec(2,:);
        phi_vec     = state_vec(3,:);
        V_vec       = state_vec(4,:)/1000;
        gamma_vec   = state_vec(5,:);
        psi_vec     = state_vec(6,:);
        m_vec       = state_vec(7,:);
        Q_vec       = state_vec(8,:);

        alpha_vec = rad2deg(alpha)*ones(length(t),1);

    % Process data
        flying_bool = 1;
        QDot_vec = 0;
        QInt_vec = 0;
        for i=1:length(r_vec)            
            % If below the earth's surface, zero everything out
            if norm(r_vec(i))-R/1000<=0
                r_vec(i)=R/1000;
                
                % Determine landing state
                if flying_bool==1
                   theta_landed = theta_vec(i);
                   phi_landed   = phi_vec(i);
                   psi_landed   = psi_vec(i);
                   m_landed     = m_vec(i);
                   Q_landed     = Q_vec(i);
                   
                   flying_bool  = 0;
                end
                
                % Set states at earth's surface properly
                theta_vec(i)   = theta_landed;
                phi_vec(i)     = phi_landed;
                V_vec(i)       = 0;
                gamma_vec(i)   = 0;
                psi_vec(i)     = psi_landed;
                m_vec(i)       = m_landed;
                Q_vec(i)       = Q_landed;
            end
            
            
            
            % Convert radius vector to cartesian value
            er(:,i) = e_r(theta_vec(i));
            r_xy(:,i) = r_vec(i)*er(:,i);
            
            
            % Determine cost value:
            if 1<i 
               QDot_vec(i) = (Q_vec(i)-Q_vec(i-1))/(t(i)-t(i-1));
               QInt_vec(i) = sum(QDot_vec(1:i));
            end
        end

    % Plot values
        figure(fig_h.Number)
        hold all
        circle(0,0,R/1000);
        hold all
        plot(r_xy(1,1),r_xy(2,1), "og",'HandleVisibility','off')
        plot(r_xy(1,end),r_xy(2,end), "or",'HandleVisibility','off')
        plot(r_xy(1,:),r_xy(2,:))
        title('Vehicle Approach to Earth')
        xlabel('x [km]')
        ylabel('y [km]')


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
        plot(t,alpha_vec)
        ylim([-1,181])
        title('Control Input vs. Time')
        xlabel('Time [s]')
        ylabel('Alpha [deg]')
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
        
        
        figure(fig_h.Number+2)
        subplot(2,1,1)
        hold all
        plot(t,QDot_vec)
        title('Heat Rate')
        xlabel('Time [S]')
        ylabel('Heat Rate [W/cm^2]')
        
        subplot(2,1,2)
        hold all
        plot(t,QInt_vec)
        title('Heat Load')
        xlabel('Time [S]')
        ylabel('Heat Load [J/cm^2]')
    
end

figure(fig_h.Number)   
hold all
legend('\alpha=0','\alpha=20','\alpha=40','\alpha=60','\alpha=80',...
        '\alpha=100','\alpha=120','\alpha=140','\alpha=160','\alpha=180')

figure(fig_h.Number+1)
subplot(2,3,4)
hold all
legend('\alpha=0','\alpha=20','\alpha=40','\alpha=60','\alpha=80',...
        '\alpha=100','\alpha=120','\alpha=140','\alpha=160','\alpha=180')
    
figure(fig_h.Number+2)
subplot(2,1,2)
hold all
legend('\alpha=0','\alpha=20','\alpha=40','\alpha=60','\alpha=80',...
        '\alpha=100','\alpha=120','\alpha=140','\alpha=160','\alpha=180')

    
%% Functions
% Radial unit vector in XY plane (for plotting)
function er = e_r(theta)
    er = [ cos(theta) ;...
            sin(theta)];
end % end e_r

% Velocity unit vector in XY plane (for plotting)
function ev = e_v(theta,fpa)
    phi_star = deg2rad(90) + fpa - theta;
    ev = [ cos(phi_star) ;...
            sin(phi_star)];
end % end e_v
    
