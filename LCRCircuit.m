%% Coded for the fulfilment of Master's Degree at Politecnico Di Milano
% Author:: Eshwar Bhargav Bhupanam
% Course:: Modeling and Simulation of Aerospace Systems
% Topic:: Electrical Circuit - LCR
% Year:: 2020-2021

%% Initializing....
clear vars; clc; close all

%% LCR Circuit - Physical model inputs
par.Vc0 = 1;     % Voltage drop in Volts 'V'
par.R1 = 1000;
par.R2 = 100;
par.L = 1e-3;   % Inductance in Henry 'H'
par.C = 1e-3;   % Capacitance in Farad 'F'
par.f = 5;      % Frequency in Hertz 'Hz'

%% Mathematical model

% Case 1
tspan = [0 3];
Vc0 = [par.Vc0 0];    % Initial value
[t_two,y] = electricalsystem (tspan,Vc0,par,0);

% Plotting
figure ()
plot(t_two,y(:,1)); grid on;
legend({'\({V_C}\) without voltage source'}, 'Interpreter', 'latex');
ylabel('Voltage \(V_c\) [V]', 'Interpreter', 'latex')
xlabel('Time [s]', 'Interpreter', 'latex')

% Case 2
tspan = [0 3];
Vc0 = [par.Vc0 0];    % Equilibrium point
[t_two,y] = electricalsystem (tspan,Vc0,par,1);

% Plotting
figure ()
plot(t_two,y(:,1)); grid on;
legend({'\({V_C}\) with voltage source'}, 'Interpreter', 'latex');  % ,'\(V_c\)'
ylabel('Voltage \(V_c\) [V]', 'Interpreter', 'latex')
xlabel('Time [s]', 'Interpreter', 'latex')

%% Functions

function [t,y] = electricalsystem (tspan,Vc0,par,type)
syms Vc(t) V(t) t
if type == 0
    V(t)=0;
else
    V = sin(2*pi*par.f*t)*atan(t);
end
EOM = (diff(Vc,2)*(par.L*par.C*(1-(par.R2/par.R1))))+...
    (((par.R2*par.C)-(par.L/par.R1))*diff(Vc))+(Vc)...
    == (V+((par.L/par.R1)*diff(V)));
EOM = simplify(EOM);
EOM = isolate(EOM,diff(Vc,2));
V = odeToVectorField(EOM);
LCRode = matlabFunction(V,'vars', {'t', 'Y'},'file','LCRode.m');
options = odeset;

[t,y] = ode45(@LCRode,tspan,Vc0,options);
delete LCRode.m
end
