

p.g = 9.81; % gravedad
p.I_sp = 250; % impulso específico
p.m0 = 1000; % masa inicial del misil
p.eps = 0.9; % eficiencia del motor
p.tb = 5; % tiempo de quemado del motor
max_step = 0.01; % paso máximo de integración
tspan = 0:max_step:1000000; % tiempo de simulación
p.Rt = 6371e3; % radio de la tierra (nota: 6371e3 m, no 6371e6)
p.Cd = 0; %[-]
p.Sref = pi*(0.1/2)^2;

y0 = [5; 10*pi/100; 0.1; 0];


options = odeset('MaxStep',max_step,'RelTol',1e-8,'AbsTol',1e-4, 'Events',@(t,y) event_function(t,y,p));

% --- 1. Simulación con Gravedad Constante ---
sol=ode45(@(t,y) Misil2D(t,y,p), tspan, y0, options);
t = sol.x;
v = sol.y(1,:);
phi = sol.y(2,:);
z = sol.y(3,:);
psi = sol.y(4,:);

X = (1+z/p.Rt).*sin(psi);
Z = (1+z/p.Rt).*cos(psi);

% Función de gravedad variable (usada en modelo 1)
function g = gravity(y, p)
    g = p.g * p.Rt^2 / (y + p.Rt)^2;
end

function dy_dt=Misil2D(t,y,p)
    dy_dt(1,1) = p.g*p.I_sp*p.eps/p.tb/(1-p.eps*t/p.tb)*(t<p.tb) -p.g/(1+y(3)/p.Rt)^2*cos(y(2));
    dy_dt(2,1) = (p.g/y(1)/(1+y(3)/p.Rt)^2-y(1)/(p.Rt+y(3)))*sin(y(2));
    dy_dt(3,1) = y(1)*cos(y(2));
    dy_dt(4,1) = y(1)/(p.Rt+y(3))*sin(y(2));
end

% Evento de impacto (sirve para ambas)
function [position, isterminal, direction] = event_function(t, y, p) 
    position = [t-p.tb; y(3); y(1)*cos(y(2))]; 
    isterminal = [0; 1; 0]; 
    direction = [-1; -1; -1]; 
end

%% --- Gráficas ---
% -- Gráfica 1 --
figure(1);
subplot(2,2,1);
plot(t,v);
ylabel('Velocidad [m/s]');
xlabel('Tiempo (s)');
grid
hold on

subplot(2,2,2);
plot(t,phi);
ylabel('Phi [rad]');
xlabel('Tiempo (s)');
grid
hold on

subplot(2,2,3);
plot(t,z/1000);
ylabel('Altura [km]');
xlabel('Tiempo (s)');
grid

subplot(2,2,4);
plot(t,psi*p.Rt/1000);
ylabel('Alcance [km]');
xlabel('Tiempo (s)');
grid
hold off

% -- Gráfica 2 --
figure(2);
theta = linspace(0,2*pi, 1000);
plot(cos(theta), sin(theta), 'k')
hold on;
plot(X, Z);
axis equal;
