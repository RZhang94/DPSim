clear; clc;
%% System Setup & IC
spatInc = 0.01;
tempInc = 0.0001;
x = 0:spatInc:1;
t = 0:tempInc:0.5;
lamb = 12; %Reaction term coefficent
%Side index to allow deriv calc
xghost = 0-spatInc:spatInc:1+spatInc; 

xdim = size(x,2)+2;
tdim = size(t,2);

%Stores values across time/space
u = zeros(xdim,tdim); 
%Initial conditions & Boundary cond
u(:,1) = 2*sin(xghost*pi); 
u(1,1) = u(2,1) - (u(3,1)-u(2,1)); 
u(xdim,1) = u(xdim-1,1) + (u(xdim-1,1)-u(xdim-2,1));
%Calculate IC control
u(:,1) = contr(u, 1, x, xdim, spatInc, lamb); 
%If no control
% ti = 1
% u(2,ti) = 0;
% u(xdim-1,ti) = 0;
% u(1,ti) = u(2,ti) - (u(3,ti)-u(2,ti));
% u(xdim,ti) = u(xdim-1,ti) + (u(xdim-1,ti)-u(xdim-2,ti));

%% Check derivative & BC
% uxreal = pi*cos(xghost*pi);
% uxxreal = -pi^2*sin(xghost*pi);
% hold on
% plot(xghost,u(:,1))
% plot(xghost,uxreal(:))
% plot(xghost,uxxreal(:))
% hold off

%% Evolution of System through time
for ti = 2:tdim
    %Evolution from previous time increment via RK4 function
    u(:,ti) = RK4(u, ti, xdim, tdim, spatInc, tempInc, lamb);
    %Define control input and reinforce boundary conditions
    u(:,ti) = contr(u, ti, x, xdim, spatInc, lamb); 
    
%   If uncontrolled reaction, comment out other block
%     u(2,ti) = 0;
%     u(xdim-1,ti) = 0;
%     u(1,ti) = u(2,ti) - (u(3,ti)-u(2,ti));
%     u(xdim,ti) = u(xdim-1,ti) + (u(xdim-1,ti)-u(xdim-2,ti));
end

%% Output Crap if you want leave, if not comment out
% meshplot(u, x, t, xdim)
animation(u,x,t,xdim,tdim)

%% Plot mesh of u
function [] = meshplot(u, x, t, xdim)
% %Create grid
[xg,tg] = meshgrid(x,t);
%Map system values to grid
u_g = u(2:xdim-1,:);
%Plot
mesh(t,x,u_g);
title('Diffusion-Reaction System with Backstepping Control')
xlabel ('t')
ylabel('x')
end

%% Animation plot of u
%Animation Plot
function [] = animation(u,x,t,xdim,tdim)
clf('reset')
a = 1;
while a == 1
    for i = 1:20:tdim   
       hold on
       grid on
       axis([0 1 -6 6]);
       h(a) = plot( x,u(2:xdim-1,i), 'b');
       title("Time: ", t(i));
       ylabel('x position');
       xlabel('spatial position');
       hold off
       drawnow;
       pause(.025);
       delete(h(a));
    end
    answer = questdlg('Replay?', 'Replay Dialog', 'Yes', 'No', 'No');
    switch answer
        case 'Yes'
            a = 1;
        case 'No'
            a = 0;
    end
end
end

%% Calculation of Control Input
function [u_c] = contr(u, ti, x, xdim, spatInc, lamb)
u_c = u(:,ti);
%Define bessel function for line
besself = x.*lamb./sqrt(lamb*(1-x.^2)).*besseli(1, sqrt(lamb*(1- x.^2))).*u_c(2:xdim-1).';
%Initial besself guess
besself(xdim-2) = 0;
%Guess the integration input step a few times to ensure accuracy of u(1,ti)
for a = 1:3
    controller = 0;
    for i= 2:size(besself,2)
        %Determine area of integral with trapezoids
        controller = controller + (besself(i) - besself(i-1))*spatInc/2 + besself(i-1)*spatInc;
    end
    %Sub in inputs and then recalculate
    besself(xdim-2) = controller;
end
controller = - controller;

%Substitute input load side
u_c(xdim-1) = controller;
%Redefine linear relation at edges
u_c(2) = 0;
u_c(1) = u_c(2) - (u_c(3)-u_c(2)) ;
u_c(xdim) = u_c(xdim-1) + u_c(xdim-1)-u_c(xdim-2);
end

%%  RK4 function
%Time in RK4 is not shifted forward as no contribution to time deriv.
function [u_n] = RK4(u, ti, xdim, tdim, spatInc, tempInc, lamb)
%Advancement from previous step
u_f = u(:,ti-1);
%Calculate k1
k1 = evolve(u_f, xdim, tdim, spatInc, tempInc, lamb);
%Calculate k2
k2 = evolve(u_f + tempInc*k1.'/2, xdim, tdim, spatInc, tempInc, lamb);
%Calculate k3
k3 = evolve(u_f + tempInc*k2.'/2, xdim, tdim, spatInc, tempInc, lamb);
%Calculate k4
k4 = evolve(u_f + tempInc*k3.'/2, xdim, tdim, spatInc, tempInc, lamb);
%Final composition
u_n = u_f.' + (k1+2*k2+2*k3+k4)*tempInc/6;
end

%% Dynamics
%The deriv wrt. time 
function [u_n] = evolve(u_o, xdim, tdim, spatInc, tempInc, lamb)
%Dynamics of system are encoded here
% u_o is the original data, u_n will be the resultant rate of change for
% next increment
% u__t = u__xx + lamb*u_o
%u(n+1) = u(n) + u__t*tempInc

%Storage vectors for derivatives
uxx = 1:xdim;

%Second derivative from central difference approx
for xi = 2:xdim-1
    uxx(xi) = (u_o(xi+1) -2*u_o(xi) +u_o(xi-1))/(spatInc^2);
end

%Ghost values at end of vector
uxx(2) = 0;
uxx(xdim-1) = 0;

u_n = uxx+ u_o.'*lamb;
end
