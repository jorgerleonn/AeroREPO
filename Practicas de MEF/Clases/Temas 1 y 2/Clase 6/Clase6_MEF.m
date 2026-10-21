%% Ejemplo de integración Gauss-Legendre de una función armónica
% Calcula la resultante de una función distribuida f_x
clear; close all; clc;

%% Valores iniciales y definición de la función
x1 = 0; 
x2 = 5; 
deltaL = x2 - x1;

C1 = 1; 
C2 = 2; 
C3 = 3;

% Declaración simbólica
syms x chi real
f_x = C1*sin(pi*x) + C2*sin(2*pi*x) + C3*sin(3*pi*x);

%% Integración numérica (Gauss-Legendre)
tStart1 = tic;

% Normalización: cambio de variable x -> chi en [-1, 1]
% x(chi) = x1 + (deltaL/2)*(chi + 1)
f_chi = subs(f_x, x, x1 + deltaL*(chi + 1)/2);

nint = 17; % Número de puntos de integración
[Xint, Wint] = GaussLegendre_WP(nint, -1, 1);

Je = deltaL / 2; % Jacobiano de la transformación
Fint_num = 0;

for ip = 1:nint
    wi = Wint(ip);
    ci = Xint(ip);
    fchi = double(subs(f_chi, chi, ci));
    Fint_num = Fint_num + fchi * wi * Je;
end

tEnd1 = toc(tStart1);

%% Integración analítica / simbólica
tStart2 = tic;
Fint_sym = double(int(f_x, x, [x1, x2]));
tEnd2 = toc(tStart2);

% Error relativo (%)
error_Fint = abs(Fint_sym - Fint_num) / abs(Fint_sym) * 100;

%% Visualización de resultados en consola
fprintf('Cuadratura de integración: %d puntos\n', nint);
fprintf('F integrada numéricamente: %.6f (Tiempo: %.4f s)\n', Fint_num, tEnd1);
fprintf('F integrada simbólicamente: %.6f (Tiempo: %.4f s)\n', Fint_sym, tEnd2);
fprintf('Error relativo: %.6e %%\n', error_Fint);

%% Representación gráfica
figure('Color', 'w');
Fx_num = matlabFunction(f_x); % Conversión a función numérica vectorizada
x_vals = linspace(x1, x2, 1000);
y_vals = Fx_num(x_vals);

area(x_vals, y_vals, 'FaceColor', [0.2 0.6 0.8], 'FaceAlpha', 0.5);
grid on;
xlabel('x');
ylabel('f(x)');

titulos = {
    sprintf('Cuadratura de integración: %d puntos', nint), ...
    sprintf('F_{num} = %.4f | F_{sym} = %.4f', Fint_num, Fint_sym), ...
    sprintf('Error relativo = %.3e %%', error_Fint)
};
title(titulos, 'FontSize', 11);

%% Subrutinas locales
function [chi_int, w_int] = GaussLegendre_WP(n, a, b)
    % Pesos y abscisas para cuadratura de Gauss-Legendre en [a, b]
    n = n - 1;
    n1 = n + 1; 
    n2 = n + 2;
    xu = linspace(-1, 1, n1)';
    
    % Estimación inicial
    y = cos((2*(0:n)' + 1)*pi / (2*n + 2)) + (0.27/n1)*sin(pi*xu*n/n2);
    
    L = zeros(n1, n2);
    Lp = zeros(n1, n2);
    
    y0 = 2;
    while max(abs(y - y0)) > eps
        L(:, 1) = 1;
        Lp(:, 1) = 0;
        L(:, 2) = y;
        Lp(:, 2) = 1;
        
        for k = 2:n1
            L(:, k+1) = ((2*k - 1)*y.*L(:, k) - (k - 1)*L(:, k-1)) / k;
        end
        
        Lp = n2 * (L(:, n1) - y.*L(:, n2)) ./ (1 - y.^2);
        y0 = y;
        y = y0 - L(:, n2) ./ Lp;
    end
    
    % Mapeo lineal a [a, b]
    chi_int = flipud((a*(1 - y) + b*(1 + y)) / 2);
    w_int = flipud((b - a) ./ ((1 - y.^2) .* Lp.^2) * (n2/n1)^2);
end