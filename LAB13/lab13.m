clc, clear all, clf

% solve u_t + p(x)*u_x = nu*u_xx; u(x,0) = u0(x)
p = @(x) sin(2*pi*x);               % sine wave
uinit = @(x) sin(2*pi*x);              % sine wave


% parameters
nu = 0.002;

% final time
ns = 128;           % number of steps
tf = 1.;            % final time
dt = tf / ns;       % step size
t = 0.;             % initial time

% create domain
np = 64;           % number of points (in space)
dx = 1. / (np);     % step size (in space)
x = 0:dx:(1-dx);    % domain
disp(x);
% sample initial condition
u0 = uinit(x);
disp(u0);
% initialize other arrays
u = zeros(length(u0),1);    % solution
A = zeros(np,np);           % for backward euler system
b = zeros(np,1);            % "

% plot initial
%plot(x,u0,'b');
%axis([0 1 -1.05 1.05]);

% time loop
for k = 1:ns

    % pause statement for controlling updates
    %pause;

    % calculate u^(n+1)
    for i = 1:np
        %------- forward euler -------%

        %%%% find my left/right neighbor (periodic)
        %%%if i==1
        %%%    ul = u0(end);
        %%%    ur = u0(i+1);
        %%%elseif i==np
        %%%    ul = u0(i-1);
        %%%    ur = u0(1);
        %%%else
        %%%    ul = u0(i-1);
        %%%    ur = u0(i+1);
        %%%end

        %%%% update point in time
        %%%unew(i) = u0(i) + dt * (-p(x(i)) * (ur - ul) / (2.*dx) ...
        %%%            + nu*(ur - 2.*u0(i) + ul) / dx^2);

        %------- backward euler -------%

        % add to matrix
        if i == 1
            A(i,i+1) = -dt * (-0.5*p(x(i))/dx + nu/dx^2) ...
                        / (1 + 2*nu*dt/dx^2);
            A(i,i) = 1.;
            A(i,end) = -dt * (0.5*p(x(i))/dx + nu/dx^2) ...
                        / (1 + 2*nu*dt/dx^2);
        elseif i == np
            A(i,1) = -dt * (-0.5*p(x(i))/dx + nu/dx^2) ...
                        / (1 + 2*nu*dt/dx^2);
            A(i,i) = 1.;
            A(i,i-1) = -dt * (0.5*p(x(i))/dx + nu/dx^2) ...
                        / (1 + 2*nu*dt/dx^2);
        else
            A(i,i+1) = -dt * (-0.5*p(x(i))/dx + nu/dx^2) ...
                        / (1 + 2*nu*dt/dx^2);
            A(i,i) = 1.;
            A(i,i-1) = -dt * (0.5*p(x(i))/dx + nu/dx^2) ...
                        / (1 + 2*nu*dt/dx^2);
        end

        % rhs
        b(i) = u0(i) / (1. + 2.*nu*dt/dx^2);

    end

    % solve system (backward euler)
    u = A \ b;

    % save new solution
    u0 = u;

    % caclulate new time
    t = t + dt;

    % plot
    %plot(x,u,'b');
    %axis([0 1 -1.05 1.05]);

end
disp(u);
disp(t);
disp(x);