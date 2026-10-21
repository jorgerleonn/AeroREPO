%% Clase 1 VA1 - Análisis Estructural por Elementos Finitos (Cerchas 2D)
clear; clc; close all;

%% 1. Parámetros y Materiales
L = 1;              % Longitud característica de la estructura [m]
E = 1;              % Módulo de elasticidad del material [Pa]
A = 1;              % Área de la sección transversal de las barras [m^2]
P = 0.025;          % Carga externa aplicada [N]
gamma = pi/4;       % Ángulo de aplicación de la carga [rad]

% Vector de fuerzas nodales globales (aplicadas en el nodo 4, grados de libertad 7 y 8)
FL = [ P*cos(gamma); 
      -P*sin(gamma)];

%% 2. Geometría y Nodos
% Coordenadas de los nodos [X, Y]
n1 = [0, L/2];
n2 = [0, 0];
n3 = [0, -L];
n4 = [L, 0];

c_nod = [n1; n2; n3; n4];           % Matriz de coordenadas globales
conectividad_n = [1 4; 2 4; 3 4];   % Definición de barras: [nodo_inicial, nodo_final]
n_elementos = size(conectividad_n, 1);
n_grados_libertad = 2 * size(c_nod, 1);

%% 3. Matrices Elementales
% Matriz de rigidez local estándar para un elemento tipo barra (cercha 2D)
M_elemento = [ 1  0 -1  0;
               0  0  0  0;
              -1  0  1  0;
               0  0  0  0];

% Inicialización de matrices globales y variables de almacenamiento
K = zeros(n_grados_libertad);   % Matriz de rigidez global
T_matrices = cell(n_elementos, 1); % Para almacenar las matrices de rotación de cada barra
L_barras = zeros(n_elementos, 1);  % Para almacenar las longitudes de cada barra

%% 4. Ensamblaje de la Matriz de Rigidez Global
for i = 1:n_elementos
    na = conectividad_n(i, 1);
    nb = conectividad_n(i, 2);
    nodo_a = c_nod(na, :);
    nodo_b = c_nod(nb, :);
    
    % Longitud de la barra
    Lei = norm(nodo_b - nodo_a);
    L_barras(i) = Lei; 
    
    % Cosenos directores
    c = (nodo_b(1) - nodo_a(1)) / Lei;
    s = (nodo_b(2) - nodo_a(2)) / Lei;
    
    % Matriz de rotación nodal y elemental
    Re = [ c, s; 
          -s, c];
    Te = [Re, zeros(2); 
          zeros(2), Re];
    
    % Guardamos la matriz de rotación para usarla luego en las tensiones
    T_matrices{i} = Te;
    
    % Rigidez elemental en ejes locales y globales
    Ke_local = ((A * E) / Lei) * M_elemento;
    Ke_global = Te' * Ke_local * Te;
    
    % Grados de libertad globales de la barra actual
    gdl_i = [2*na-1, 2*na, 2*nb-1, 2*nb];
    
    % Ensamblaje en la matriz global
    K(gdl_i, gdl_i) = K(gdl_i, gdl_i) + Ke_global; 
end

%% 5. Resolución del Sistema (Desplazamientos)
% Identificación de grados de libertad (Nodos 1, 2, 3 están fijos; Nodo 4 es libre)
gdl_libres = [7, 8]; 
K_LL = K(gdl_libres, gdl_libres);

% Cálculo de desplazamientos (F = K*u -> u = K\F)
u_l = K_LL \ FL; 

% Reconstrucción del vector de desplazamientos global completo
u_global = zeros(n_grados_libertad, 1);
u_global(gdl_libres) = u_l;

% Cálculo de las nuevas coordenadas (deformadas)
c_nod_def = c_nod;
for i = 1:size(c_nod, 1)
    c_nod_def(i, :) = c_nod(i, :) + [u_global(2*i - 1), u_global(2*i)];
end

%% 6. Cálculo de Tensiones
tensiones_barras = zeros(n_elementos, 1);

for i = 1:n_elementos
    na = conectividad_n(i, 1);
    nb = conectividad_n(i, 2);
    gdl_i = [2*na-1, 2*na, 2*nb-1, 2*nb]; 
    
    % Extraemos desplazamientos globales de la barra y pasamos a locales
    u_barra_global = u_global(gdl_i);
    u_barra_local = T_matrices{i} * u_barra_global;
    
    % Desplazamientos axiales (locales) en los nodos A y B
    u1_local = u_barra_local(1); 
    u2_local = u_barra_local(3); 
    
    % Deformación y tensión (Ley de Hooke)
    deformacion = (u2_local - u1_local) / L_barras(i);
    tensiones_barras(i) = E * deformacion;
end

% Mostrar resultados por consola
disp('--- RESULTADOS ---');
disp('Desplazamientos del nodo 4 [ux, uy]:');
disp(u_l');
disp('Tensiones en cada barra (Positivo = Tracción, Negativo = Compresión):');
for i = 1:n_elementos
    fprintf('Barra %d: %8.4f Pa\n', i, tensiones_barras(i));
end

%% 7. Representación Gráfica
figure('Name', 'Estructura Original vs Deformada', 'Color', 'w');
hold on; grid on; axis equal;

% Dibujar estructura original (Negro)
for i = 1:n_elementos
    na = conectividad_n(i, 1);
    nb = conectividad_n(i, 2);
    plot([c_nod(na, 1), c_nod(nb, 1)], [c_nod(na, 2), c_nod(nb, 2)], 'k-o', 'LineWidth', 1.5);
end

% Dibujar estructura deformada (Rojo)
for i = 1:n_elementos
    na = conectividad_n(i, 1);
    nb = conectividad_n(i, 2);
    plot([c_nod_def(na, 1), c_nod_def(nb, 1)], [c_nod_def(na, 2), c_nod_def(nb, 2)], 'r-*', 'LineWidth', 1.5);
end

title('Clase 1 y 2 - Deformación de la Estructura');
xlabel('Posición X [m]');
ylabel('Posición Y [m]');
legend('Original', 'Deformada', 'Location', 'best');
hold off;