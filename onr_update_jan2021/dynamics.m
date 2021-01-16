%Skye Mceowen
%3D Projectile Dynamics - ONR Update Sims
%Jan15, 2021

function xDot = dynamics(t,x,sigma,alpha,T)

    r       = x(1);
    theta   = x(2);
    phi     = x(3);
    V       = x(4);
    gamma   = x(5);
    psi     = x(6);
    m       = x(7);
    
    % Constants
    g0 = 9.81; %m/s^2
    G = 6.67408e-11; %[m^3/(kg s^2)], gravitational constant
    M = 5.972e24; % [kg], mass of the earth
    R = 6371e3; %[m], radius of the earth
    mu = G*M; %[m], gravitational standard param
    Isp = 250; % [s], specific impulse
    d = 1.4; % [m], body diameter
    A = (pi*(d^2))/4; %[m^2], cross-sectional area of waverider
    %m = 1200; %[kg]
    
    
    % Rho, source: https://www.spaceacademy.net.au/watch/debris/atmosmod.htm
        % TODO: use ONR model-->1976 ICAO model
    rho0 = 1.3; %[kg/m^2], 
    H = 7000; 
    h = norm(r) - R;
    if norm(r) >= R
        rho = rho0*exp(-h/H); %kg/m^3, density; ISOTHERMAL MODEL
    else
        rho = rho0;
    end
    
    % Gravity
    g = mu/(r^2);

    % Lift and drag ratio coefficients
    Cl = -0.04 + 0.8*alpha;
    Cd = 0.012 - 0.01*alpha + 0.6*alpha^2;
    
    % Lift 
    L = 0.5*Cl*rho*(V^2)*A;
    
    % Drag
    D = 0.5*Cd*rho*V^2*A;

    % Determine derivatives 
    rDot = V*sin(gamma);
    
    thetaDot = (V*cos(gamma)*sin(psi))/(r*cos(phi));
    
    phiDot = (V*cos(gamma)*cos(psi))/r;
    
    vDot = (T/m)*cos(alpha) - (D/m) - g*sin(gamma);
    
    gammaDot = ( ( (T/m)*sin(alpha) + L/m )*cos(sigma) - (g - ((V^2)/r))*cos(gamma) )/V;
    
    psiDot = ( ( (T/m)*sin(alpha)  + L/m )*(sin(sigma)/cos(gamma)) + (V^2*cos(gamma)*sin(psi)*tan(phi))/r )/V;
    
    mDot = -T/(Isp*g0);
    
    % Stack into state propogation vector
    xDot = [rDot; thetaDot; phiDot; vDot; gammaDot; psiDot; mDot];

end