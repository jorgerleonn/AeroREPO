%% Introducicción a MATLAB
% B = [1 2 3; 4 5 6; 7 8 9];
% disp(B);

% representación de la función seno de forma discreta
% v = linspace(0, 2*pi, 100);
% disp(v);
% y = sin(v);
% plot(v, y);

%% ---------------------- Código ode45 Clase 1 - Movimiento 1D de un misil --------------------

params.g = 9.81; % gravedad
params.I_sp = 250; % impulso específico
params.m0 = 1000; % masa inicial del misil
params.eps = 0.9; % eficiencia del motor
params.tb = 5; % tiempo de quemado del motor
max_step = 0.1; % paso máximo de integración
tspan = 0:max_step:10000; % tiempo de simulación
params.Rt = 6371e6; % radio de la tierra

options = odeset('MaxStep',max_step,'RelTol',1e-4,'AbsTol',1e-4, 'Events',@(t,y) event_function(t,y,params));

m = ode45(@(t,y)misil_1D(t,y,params), [0 tspan], [0; 0], options);

function g = gravity(y, params)
    g = params.g * params.Rt^2 / (y + params.Rt)^2;
end

function dY_dt = misil_1D(t, Y, p) % lanzamiento vertical y sin resistencia del aire
    dY_dt(1,1) = p.g*p.I_sp*p.eps/p.tb/(1-p.eps*t/p.tb)*(t<p.tb) - gravity(y(2), p); % velocidad
    dY_dt(2,1) = Y(1); % posición
end

function [position, isterminal, direction] = event_function(t, y, params) % evento para detectar cuando el misil toca el suelo
    position = [y(2), params.tb-t, y(1)]; % evento cuando la posición es cero (misil toca el suelo)
    isterminal = [1, 0, 0]; % detener la integración
    direction = [-1, -1, -1]; % solo detectar cuando la posición disminuye
end

t = m.x; % tiempo de simulación
y = m.y; % resultados de la simulación
v = y(1,:); % posición
z = y(2,:); % velocidad
vb = m.ye(1,1); % velocidad al final del quemado

figure(1);
subplot(2,1,1);
plot(t, v);
% punto de fin de combustión
plot(t, vb, 'ro', 'MarkerSize', 2, 'LineWidth', 2); % marcador rojo para la velocidad al final del quemado
ylabel('Velocidad (m/s)');
xlabel('Tiempo (s)');
title('Velocidad del misil en función del tiempo');


figure(2);
subplot(2,1,1);
plot(t, z);
ylabel('Posición (m)');
xlabel('Tiempo (s)');
title('Posición del misil en función del tiempo');