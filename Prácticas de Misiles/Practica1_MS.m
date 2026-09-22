%% Introducicción a MATLAB
% B = [1 2 3; 4 5 6; 7 8 9];
% disp(B);

% representación de la función seno de forma discreta
% v = linspace(0, 2*pi, 100);
% disp(v);
% y = sin(v);
% plot(v, y);

%% ---------------------- Código ode45 Clase --------------------

params.g = 9.81; % gravedad
params.I_sp = 250; % impulso específico
params.m0 = 1000; % masa inicial del misil
params.eps = 0.9; % eficiencia del motor
params.tb = 25; % tiempo de quemado del motor
max_step = 0.1; % paso máximo de integración
tspan = 0:max_step:50; % tiempo de simulación

options = odeset('MaxStep',max_step,'RelTol',1e-4,'AbsTol',1e-4, 'Events',@(t,y) event_function(t,y,params));

m = ode45(@(t,y) misil_1D(t,y,params), [0 tspan], [0; 0], options);

function dY_dt = misil_1D(t, Y, params) % lanzamiento vertical y sin resistencia del aire
    dY_dt(1,1) = params.g*params.I_sp*params.eps/params.tb/(1-params.eps*t/params.tb)*(t<params.tb) - params.g; % velocidad
    dY_dt(2,1) = Y(1); % posición
end

function [position, isterminal, direction] = event_function(t, y, params)
    position = [y(2)]; % evento cuando la posición es cero (misil toca el suelo)
    isterminal = [1; 0]; % detener la integración
    direction = [-1; 0]; % solo detectar cuando la posición disminuye
end

t = m.x; % tiempo de simulación
y = m.y; % resultados de la simulación
z = y(1,:); % posición
v = y(2,:); % velocidad