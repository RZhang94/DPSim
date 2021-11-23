clear; clc;
%% Physical Parameters of System
spatInc = 1/25;
tempInc = 0.0001;

%% Physical Paramters of String
rho = 0.05;
lmax = 1;
m_p = 0.25; %Mass at end of string
g = 9.81;

%% %% Physical Parameters of Drone
m = 1.6; %drone weight in kg
Jx = 0.03; %Moment of inertia about the planar axis

%% Setup
%Parameters
x = 0:spatInc:lmax;
t = 0:tempInc:20;
xghost = 0-spatInc:spatInc:lmax+spatInc;
xdim = size(x,2)+2;
tdim = size(t,2);

%Arrays
u = zeros(xdim,tdim); %x vector
ut = zeros(xdim,tdim); %x vel vector
pz = zeros(1,tdim); %z pos
vz = zeros(1,tdim); %z velocity
theta = zeros(1,tdim); %Angle in radians
omega = zeros(1,tdim); %Angular velocity

%% Initial Conditions
%Pend IC
u(:,1) = 1.2 + (xghost-lmax)*0.25; 
ut(:,1) = 1.2 + (xghost-lmax)*0.25;
%Drone IC
pz(1) = 1;
vz(1) = 1;
theta(1) = 0.125;
omega(1) = 0.125;

%% Evolution /w Time via RK4
for ti = 2:tdim
%Evolution of time increment
[u(:,ti), ut(:,ti), pz(1,ti), vz(1,ti), theta(1,ti), omega(1,ti)] =  RK4(u(:,ti-1), ut(:,ti-1),pz(1,ti-1),vz(1,ti-1),theta(1,ti-1), omega(1,ti-1), xdim,  spatInc, tempInc, xghost, m_p, rho, m, Jx);
end
%% Display Functions
animatechoice(u,xdim,pz,theta,lmax,tdim,x,t);

function [] = animatechoice(u,xdim,pz,theta,lmax,tdim,x,t)
    a = 1;
    while a == 1
    answer = questdlg('Animation Type:' , 'Choice', 'PDE Mesh Plot', 'Animation', 'Error Plots', 'Animation');
    switch answer
        case 'PDE Mesh Plot'
            meshplot(u,t,x);
        case 'Animation'
            spd = inputdlg('Frames per refresh','Input',[1,40],{'700'});
            spd = floor(str2double(spd));
            animateplot(u, pz,theta, lmax, xdim,tdim,x,t,spd);
        case 'Error Plots'
            posError(u, pz, xdim, tdim, t);
    end
    answer = questdlg('Reselect?: ', '', 'Yes', 'No', 'No');
    switch answer
        case 'Yes'
            a = 1;
        case 'No'
            a = 0;
    end
    end
end

%Plot out a surface for the evolution of the PDE system
function [] = meshplot(u,t ,x)
% mesh(u(2:xdim-1,1:50:end))
[tgr, xgr] = meshgrid(t,x);
tgr = tgr(:, 1:10:end);
xgr = xgr(:, 1:10:end);
ugr = u(2:end-1, 1:10:end);
mesh(tgr, xgr, ugr)
title('Mesh Plot of Rope Position')
xlabel('Time (t)')
ylabel('Spatial Position (s)')
zlabel('X Position')
end

%Animate what the drone simulation looks like
function [] = animateplot(u, pz,theta, lmax, xdim,tdim,x,t,spd)
clf('reset');
arm = 0.25;
a = 1;
while a == 1;
    st_t = clock;
    for i = 1:spd:tdim          
       hold on
       plot ([0 0 -50 50], [-50 50 0 0], 'r--');
       h(a) = plot(u(2:xdim-1,i), pz(i)+(x-lmax) , 'b'); 
       g(a) = plot(u(xdim-1,i), pz(i), 'ro'); 
       j(a) = plot([u(xdim-1,i)-sin(theta(i)+pi/2)*arm, u(xdim-1,i)+sin(theta(i)+pi/2)*arm], [ pz(i)-cos(theta(i)+pi/2)*arm, pz(i)+cos(theta(i)+pi/2)*arm], 'g');
       k(a) = plot(u(2, i), pz(i)-lmax, 'go');
       mid = u(floor(xdim/2)+1,i);
       axis( [u(xdim-1,i)-1 u(xdim-1,i)+1 -1.5+pz(i) 0.5+pz(i)]);
       title("Time: ", t(i));
       xlabel('x position');
       ylabel('spatial position');
       hold off
       drawnow;
       pause(.02);
       delete(h(a));
       delete(g(a));
       delete(j(a));
       delete(k(a));
    end
    fi_t = clock;
    dur = fi_t-st_t;
    answer = questdlg(strcat('Replay? Duration: ', string(dur(6))) , 'Replay Dialog', 'Yes', 'No', 'No');
    switch answer
        case 'Yes'
            a = 1;
        case 'No'
            a = 0;
    end
end
end

%Placeholder for x and z positional error function
function [] = posError(u, pz, xdim, tdim, t);
 xerror = 1:tdim;
 xerror = xerror .*0;
    for a = 1:tdim
        for b = 2:xdim-1
            xerror(a) = xerror(a) + (u(b,a)^2)^(1/2);
        end
    end
 xintError = xerror .*0;
 %Interpolate as square - Performance metric not required for calc
 spatInc = t(2)-t(1);
for a = 2:tdim
    xintError(a) = xintError(a-1) + xerror(a-1)*spatInc;
end
disp(xintError(tdim))
    
subplot(1,2,1)
plot(t, xerror)
title('X Position Error')
xlabel('Time')
ylabel('Error')

subplot(1,2,2)
plot(t, xintError)
title('Cumulative X Position Error')
xlabel('Time')
ylabel('Cumulative Error')
end


%% Functions for Deriv

function [uLB] = LBound(u_ti, xdim, spatInc)
%Input: Magnitude of DPS @ time increm ent
%Output: DPS with extra values past boundary, linearized to allow der. aprox.
uLB = u_ti;
uLB(1) = uLB(2) - (uLB(3)-uLB(2)) ;
uLB(xdim) = uLB(xdim-1) +( uLB(xdim-1)- uLB(xdim-2));
end

function [ux, uxx] = Derivs(u_ti, xdim, spatInc)
%Input: Magnitude of DPS in vector @ time increment
%Output: Vectors for first and second derivative approximation
u = LBound(u_ti,xdim,spatInc);
ux = u.*0;
uxx = u.*0;
for xi = 2:xdim-1
   ux(xi) = ( u(xi+1) - u(xi-1)) / (2*spatInc);
end
ux = LBound(ux, xdim, spatInc);
%Second deriv wrt x
for xi = 2:xdim-1
   uxx(xi) = (u(xi+1) - 2*u(xi) + u(xi-1)) / (spatInc)^2;
end
end

%% RK4 function
%Time in RK4 is not shifted forward as no contribution to time deriv.
%Both RK4s need to run in parallel
function [u_n, ut_n, pz_n, vz_n, theta_n, omega_n] = RK4(u_ti, ut_ti, pz, vz, theta, omega, xdim, spatInc, tempInc, xghost, m_p, rho, m, Jx)
%Advancement from previous step

%Construct duples for resources
xghost = xghost.';
u = [LBound(u_ti,xdim,spatInc), LBound(ut_ti,xdim,spatInc)];
d = [pz, vz, theta, omega];
p = [xdim, spatInc, tempInc, m_p, rho, m, Jx, 9.81];

%Calculate k1:k4
[k1a, k1b, k1c, k1d, k1e, k1f, u2, d2] = DynFunc(u,p,d,xghost, 2);
[k2a, k2b, k2c, k2d, k2e, k2f, u3, d3] = DynFunc(u2,p,d2,xghost, 2);
[k3a, k3b, k3c, k3d, k3e, k3f, u4, d4] = DynFunc(u3,p,d3,xghost, 1);
[k4a, k4b, k4c, k4d, k4e, k4f, u5, d5] = DynFunc(u4,p,d4,xghost, 0);

%Final composition, temp included already
u_n = u_ti + (k1a+2*k2a+2*k3a+k4a)/6;
pz_n = pz + (k1b+2*k2b+2*k3b + k4b)/6;
ut_n = ut_ti + (k1c+2*k2c+2*k3c+k4c)/6;
vz_n = vz +(k1d+2*k2d+2*k3d + k4d)/6;
theta_n = theta+ (k1e+2*k2e+2*k3e + k4e)/6;
omega_n = omega+ (k1f+2*k2f+2*k3f + k4f)/6;

end

%% RK4 Components
function [kna, knb, knc, knd, kne, knf, u_n, d_n] = DynFunc(u,p,d,xghost, fac)
%Calculate new dynamics for time period requested. Factor should be 1,2,2,1
%Calculate controller inputs for this iteration in RK4
F = 0;
%Calculate increment change based on system dynamics
kna = fx1(u,p)*p(3);
knb = 0;
knc = fx3(u,d,p,xghost,F)*p(3);
knd = 0;
kne = 0;
knf = 0;

%Calculate resultant vectors for next RK4 iteration
u_n = u + [kna, knc]./fac;
u_n = [LBound(u_n(:,1), p(1), p(2)), LBound(u_n(:,2), p(1),p(2))];
d_n = d + [knb, knd, kne, knf]./fac;
end

%% Dynamics
%Dynamics of x position
function [fx1] = fx1(u,p)
%Input: Previous increments velocity Vector
%Output: Change in x position Vector (vel)
fx1 = u(:,2);
fx1(p(1)-1) = RefGen(u,p);
end
%Dynamics of z position
function [fx2] = fx2(d)
%Input: Previous Increments Velocity in z position
%Output: Change in z position (vel)
fx2 = d(2);
end
%Dynamics of x velocity
function [fx3] = fx3(u,d,p,xghost, F)
%Input: Previous increments x Vectors
%Output: Change in x velocity Vector 
[ux, uxx] = Derivs(u(:,1), p(1), p(2)); %calculate derivs
Fx = F*sin(d(3));

fx3 = p(8).*ux + p(8).*xghost.*uxx + p(4)*p(8).*uxx./p(5);     
fx3(p(1)-1) = 0; %Velocity controlled
fx3(2) = p(8)*ux(2);
end
%Dynamics of z velocity
function [fx4] = fx4(u,d,p, F)
%Input: Previous increments z vel
%Output: change in z vel
Fz = F*cos(d(3));
sy_m = (p(1)-3)*p(2)*p(5)+ p(4) + p(6);
fx4 = Fz/sy_m - p(8);
end
%Dynamics of drone rotation
function [fx5] = fx5(d)
%Input: Previous increments angle accel
%Output: Angular velocity
fx5 = d(4);
end
%Dynamics of drone rotational velocity
function [fx6] = fx6(u, d, p, Tor)
%Input: Previous increments angular vel
%Output: Angular accel
fx6 = Tor/p(7);
end

%% Controllers
%Sub function to generate reference velocity, Kp/Kd gain in here
function [vel_ref] = RefGen(u, p)
Kp = 1.5;  %Gain matrix of Lyapunov based feedback
Kd = 1 ;
[ux, uxx] = Derivs(u(:,1), p(1), p(2));
vel_ref = -Kd * ux(p(1)-1)- Kp* u(p(1)-1,1)+((u(2,2)*ux(2)*Kd*p(4)*p(8)/p(5)/(Kd * ux(p(1)-1) * p(8) * ( 1 + p(8)/p(5))+ Kp* u(p(1)-1,1))));
if isnan(vel_ref)
    vel_ref = 0;
end
end

%Generate control inputs based on error of reference velocity
function [Tor, F] = InGen(u, p , d)
vel_ref = RefGen(u,p); %Call sub

%Gain matrix for drone seeking velocity
K = [0,1.5,0,3,0,0; 0,0,0.5,0,1.5,0.25]; 

errArr = [u(p(1)-1,1)-0, d(1)-0, u(p(1)-1,2)-vel_ref, d(2)-0, d(3)-0, d(4)-0];
F = -K(1,:)*errArr.' + ((p(1)-3)*p(2)*p(5)+ p(4) + p(6))*p(8);
Tor = -K(2,:)*errArr.';
end