%% Placa Rayleigh-Ritz con 2 funciones ortogonales (x, y)
% Teoría de flexión de placas de Mecánica de Sólidos
% Material ortótropo en direcciones fuera de la simetría de ortotropía (anisótropo aparente)
clear; close all; clc;

%% Geometría y variables simbólicas
a = 1000;      % Longitud en x [mm]
b = a / 2;     % Longitud en y [mm]
t = a / 100;   % Espesor [mm]

syms x y z chi eta real

%% Propiedades del material ortótropo 2D
E1 = 150e3;    % [MPa]
E2 = 10e3;     % [MPa]
nu12 = 0.35;
G12 = 5e3;     % [MPa]

nu21 = nu12 * E2 / E1;

% Rigideces a flexión en ejes principales de ortotropía
D11 = (t^3 / 12) * E1 / (1 - nu12 * nu21);
D22 = (t^3 / 12) * E2 / (1 - nu12 * nu21);
D12 = nu12 * D22;
D66 = (t^3 / 12) * G12;
Dtp_orto = [D11, D12, 0; 
            D12, D22, 0; 
            0,   0,   D66];

%% Rotación fuera de los ejes de simetría (Giro theta = 30º)
theta = 30 * pi / 180;
c = cos(theta); s = sin(theta);
c2 = c^2; s2 = s^2;
c3 = c^3; s3 = s^3;
c4 = c^4; s4 = s^4;
s2c2 = s2 * c2;

% Matriz de rigideces transformada (comportamiento anisótropo acoplado)
D11a = D11*c4 + 2*(D12 + 2*D66)*s2c2 + D22*s4;
D12a = (D11 + D22 - 4*D66)*s2c2 + D12*(s4 + c4);
D22a = D11*s4 + 2*(D12 + 2*D66)*s2c2 + D22*c4;
D16a = (D11 - D12 - 2*D66)*(s*c3) + (D12 - D22 + 2*D66)*(s3*c);
D26a = (D11 - D12 - 2*D66)*(s3*c) + (D12 - D22 + 2*D66)*(s*c3);
D66a = (D11 + D22 - 2*D12 - 2*D66)*s2c2 + D66*(s4 + c4);

Dtp_ani = [D11a, D12a, D16a;
           D12a, D22a, D26a;
           D16a, D26a, D66a];

%% Cargas transversales en Z
f0z = -1;                          % Presión máxima [N/mm^2]
fz_xy = f0z * (x / a) * (y / b);   % Presión bilineal
Fz = double(int(int(fz_xy, x, 0, a), y, 0, b));

%% Funciones de forma de Rayleigh-Ritz
m = 3; % Grados de aproximación por dirección (m x m funciones)
Cx = sym('Cx', [m, 1]); 
Cy = sym('Cy', [m, 1]);
Nx = sym('Nx', [1, m]); 
Ny = sym('Ny', [1, m]);

for i = 1:m
    Nx(i) = sin(pi * x * i / a); % Simplemente apoyada en x
    Ny(i) = sin(pi * y * i / b); % Simplemente apoyada en y
end

N00 = simplify((Nx * Cx) * (Ny * Cy));

NRR = sym(zeros(1, m * m));
k = 0;
for i = 1:m
    for j = 1:m
        k = k + 1;
        NRR(1, k) = diff(diff(N00, Cy(j)), Cx(i));
    end
end
dim = length(NRR);

%% Operadores diferenciales y matriz de rigidez no integrada
Nd2x  = diff(NRR, x, 2); 
Nd2y  = diff(NRR, y, 2); 
Nd2xy = diff(NRR, x, y);

% Matriz de rigidez elemental (energía de deformación por flexión)
Nd2sint = (Nd2x' * Nd2x) * D11a + ...
          (Nd2y' * Nd2y) * D22a + ...
          2 * D12a * (Nd2x' * Nd2y) + ...
          4 * D66a * (Nd2xy' * Nd2xy) + ...
          2 * D16a * ((Nd2x' * Nd2xy) + (Nd2xy' * Nd2x)) + ...
          2 * D26a * ((Nd2y' * Nd2xy) + (Nd2xy' * Nd2y));

% Vector de fuerzas nodales equivalentes sin integrar
Fzeqsint = NRR' * fz_xy;

%% Integración numérica 2D (Gauss-Legendre)
tStart1 = tic;
K0_num = zeros(dim); 
Fzeq_num = zeros(dim, 1);

nint = 10; % Cuadratura de Gauss por dimensión (nint x nint puntos)

% Cambio al dominio canónico [-1, 1] x [-1, 1]
Fzeqsint_n = subs(Fzeqsint, [x, y], [a * (chi + 1) / 2, b * (eta + 1) / 2]);
Nd2sint_n  = subs(Nd2sint,  [x, y], [a * (chi + 1) / 2, b * (eta + 1) / 2]);

[Xint, Wint] = GaussLegendre_2D(nint);
Je = (a / 2) * (b / 2); % Jacobiano 2D constante

for ip = 1:(nint * nint)
    w_i = Wint(ip, :);
    chi_ip = Xint(ip, 1); 
    eta_ip = Xint(ip, 2);
    
    weight = w_i(1) * w_i(2) * Je;
    Nd2_eval = double(subs(Nd2sint_n, [chi, eta], [chi_ip, eta_ip]));
    Fzeq_eval = double(subs(Fzeqsint_n, [chi, eta], [chi_ip, eta_ip]));
    
    K0_num   = K0_num   + Nd2_eval * weight;
    Fzeq_num = Fzeq_num + Fzeq_eval * weight;
end
tEnd1 = toc(tStart1);

%% Integración simbólica exacta
tStart2 = tic;
Fzeq_sym = double(int(int(Fzeqsint, x, 0, a), y, 0, b));
K0_sym   = double(int(int(Nd2sint, y, 0, b), x, 0, a));
tEnd2 = toc(tStart2);

%% Cálculo y presentación de errores
errk = max(abs(K0_sym(:) - K0_num(:))) / max(abs(K0_sym(:))) * 100;
errf = max(abs(Fzeq_sym - Fzeq_num)) / max(abs(Fzeq_sym)) * 100;

fprintf('--- RESULTADOS INTEGRACIÓN ---\n');
fprintf('Tiempo integración numérica: %.4f s\n', tEnd1);
fprintf('Tiempo integración simbólica: %.4f s\n', tEnd2);
fprintf('Error relativo matriz K: %.4e %%\n', errk);
fprintf('Error relativo vector F: %.4e %%\n', errf);

%% Resolución del sistema lineal y desplazamientos
C_1 = K0_num \ Fzeq_num; % Coeficientes Rayleigh-Ritz
u_z_1 = NRR * C_1;

uz1_center = double(subs(u_z_1, [x, y], [a / 2, b / 2]));
fprintf('Desplazamiento en el centro (a/2, b/2): %.4f mm\n', uz1_center);

%% Representación gráfica 3D
n_pts = 40;
x_vec = linspace(0, a, n_pts);
y_vec = linspace(0, b, n_pts);
[X, Y] = meshgrid(x_vec, y_vec);

uz_func = matlabFunction(u_z_1, 'Vars', [x, y]);
uz_eval = uz_func(X, Y);

peak = max(abs(uz_eval(:)));
escala = (0.25 * a) / peak;

figure('Color', 'w');
surf(X, Y, escala * uz_eval, uz_eval, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
colormap('turbo');
colorbar;
grid on;
axis equal;
view(-160, 26);

xlabel('X [mm]');
ylabel('Y [mm]');
zlabel('Desplazamiento escalado [mm]');
title('Deformada elástica de la placa (Rayleigh-Ritz)', 'FontSize', 12);
subtitle(sprintf('w_{max}: %.4f mm | Escala visual: %.1fx', peak, escala), 'FontSize', 10);

%% Subrutinas locales
function [P_int, W_int] = GaussLegendre_2D(nint)
    [Xint, Wint] = GaussLegendre_WP(nint, -1, 1);
    n_int = nint * nint;
    P_int = zeros(n_int, 2);
    W_int = zeros(n_int, 2);
    
    [gx, gy] = meshgrid(Xint, Xint);
    [wx, wy] = meshgrid(Wint, Wint);
    
    P_int(:, 1) = gx(:);
    P_int(:, 2) = gy(:);
    W_int(:, 1) = wx(:);
    W_int(:, 2) = wy(:);
end

function [chi_int, w_int] = GaussLegendre_WP(n, a, b)
    n = n - 1;
    n1 = n + 1; 
    n2 = n + 2;
    xu = linspace(-1, 1, n1)';
    
    y = cos((2*(0:n)' + 1)*pi / (2*n + 2)) + (0.27/n1)*sin(pi*xu*n/n2);
    
    L  = zeros(n1, n2);
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
    
    chi_int = flipud((a*(1 - y) + b*(1 + y)) / 2);
    w_int   = flipud((b - a) ./ ((1 - y.^2) .* Lp.^2) * (n2 / n1)^2);
end