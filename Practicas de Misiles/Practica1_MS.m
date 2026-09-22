%% ---------------------- Código ode45 - Movimiento 1D de un misil --------------------

params.g = 9.81; % gravedad
params.I_sp = 250; % impulso específico
params.m0 = 1000; % masa inicial del misil
params.eps = 0.9; % eficiencia del motor
params.tb = 5; % tiempo de quemado del motor
max_step = 0.1; % paso máximo de integración
tspan = 0:max_step:10000; % tiempo de simulación
params.Rt = 6371e3; % radio de la tierra (nota: 6371e3 m, no 6371e6)

options = odeset('MaxStep',max_step,'RelTol',1e-4,'AbsTol',1e-4, 'Events',@(t,y) event_function(t,y,params));

% --- 1. Simulación con Gravedad Variable ---
m_var = ode45(@(t,y)misil_1D_gvar(t,y,params), tspan, [0; 0], options);
t_var = m_var.x; 
v_var = m_var.y(1,:); 
z_var = m_var.y(2,:); 

% --- 2. Simulación con Gravedad Constante ---
m_const = ode45(@(t,y)misil_1D_gconst(t,y,params), tspan, [0; 0], options);
t_const = m_const.x; 
v_const = m_const.y(1,:); 
z_const = m_const.y(2,:); 

% Identificar el punto de fin de quemado (t = tb) para la primera simulación
[~, idx_tb] = min(abs(t_var - params.tb));
t_burn = t_var(idx_tb);
v_burn = v_var(idx_tb);

%% --- Gráficas ---

figure(1);
plot(t_var, v_var, 'b', 'LineWidth', 1.5); % Línea azul continua
hold on; 
plot(t_const, v_const, 'r--', 'LineWidth', 1.5); % Línea roja discontinua
plot(t_burn, v_burn, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6); % Punto negro
hold off;
ylabel('Velocidad (m/s)');
xlabel('Tiempo (s)');
title('Velocidad del misil en función del tiempo');
legend('Gravedad Variable', 'Gravedad Constante', 'Fin de combustión', 'Location', 'best');

figure(2);
plot(t_var, z_var, 'b', 'LineWidth', 1.5);
hold on;
plot(t_const, z_const, 'r--', 'LineWidth', 1.5);
hold off;
ylabel('Posición (m)');
xlabel('Tiempo (s)');
title('Posición del misil en función del tiempo');
legend('Gravedad Variable', 'Gravedad Constante', 'Location', 'best');


%% ---------------- Funciones Auxiliares ----------------

% Función de gravedad variable (usada en modelo 1)
function g = gravity(y, params)
    g = params.g * params.Rt^2 / (y + params.Rt)^2;
end

% Ecuaciones con gravedad variable
function dY_dt = misil_1D_gvar(t, Y, p) 
    dY_dt = zeros(2,1); 
    dY_dt(1,1) = p.g*p.I_sp*p.eps/p.tb/(1-p.eps*t/p.tb)*(t<p.tb) - gravity(Y(2), p); 
    dY_dt(2,1) = Y(1); 
end

% Ecuaciones con gravedad constante (NUEVA)
function dY_dt = misil_1D_gconst(t, Y, p) 
    dY_dt = zeros(2,1); 
    % Restamos directamente la constante p.g
    dY_dt(1,1) = p.g*p.I_sp*p.eps/p.tb/(1-p.eps*t/p.tb)*(t<p.tb) - p.g; 
    dY_dt(2,1) = Y(1); 
end

% Evento de impacto (sirve para ambas)
function [position, isterminal, direction] = event_function(t, y, params) 
    position = [y(2), params.tb-t, y(1)]; 
    isterminal = [1, 0, 0]; 
    direction = [-1, -1, -1]; 
end