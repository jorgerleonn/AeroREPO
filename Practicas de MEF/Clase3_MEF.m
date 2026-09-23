%% Clase 2 VA2 - Añadir nodos en una parte de la estructura
clear; clc; close all;

%% 1. Parámetros y Materiales
L = 1;              % Longitud característica de la estructura [m]
E = 1;              % Módulo de elasticidad del material [Pa]
A = 1;              % Área de la sección transversal de las barras [m^2]
P = 0.25;          % Carga externa aplicada [N]
gamma = -pi/4;       % Ángulo de aplicación de la carga [rad]


%% 2. Geometría y Nodos
% Coordenadas de los nodos [X, Y]
n1 = [0, L];
n2 = [0, 0];
n3 = [0, -L];
n4 = [L, 0];
n5 = [L/2, L/2];  % Nodo adicional

c_nod = [n1; n2; n3; n4; n5];           % Matriz de coordenadas globales
conectividad_n = [1 5; 2 4; 3 4; 5 4];   % Definición de barras: [nodo_inicial, nodo_final]
n_elementos = size(conectividad_n, 1); % Número de elementos barra
gdlxnod = 2; % Grados de libertad por nodo (2D)
n_grados_libertad = gdlxnod * size(c_nod, 1); % Número de grados de libertad


%% 3. Vector de fuerzas nodales globales
F = zeros(n_grados_libertad, 1);  % Inicialización del vector de fuerzas
gdl_P = [7, 8];  % Grados de libertad donde se aplica la fuerza (nodo 4)

F(gdl_P(1)) = P * cos(gamma);  % Componente X de la fuerza
F(gdl_P(2)) = P * sin(gamma);  % Componente Y de la fuerza

%% 4. Matrices Elementales
% Matriz de rigidez local estándar para un elemento tipo barra (cercha 2D)
M_elemento = [ 1  0 -1  0;
               0  0  0  0;
              -1  0  1  0;
               0  0  0  0];

% Inicialización de matrices globales y variables de almacenamiento
K = zeros(n_grados_libertad);   % Matriz de rigidez global
T_matrices = cell(n_elementos, 1); % Para almacenar las matrices de rotación de cada barra
L_barras = zeros(n_elementos, 1);  % Para almacenar las longitudes de cada barra


%% 5. Ensamblaje de la Matriz de Rigidez Global
c = zeros(n_elementos, 1); % Cosenos directores
s = zeros(n_elementos, 1); % Senos directores
Re = cell(n_elementos, 1); % Matrices de rotación nodal

for i = 1:n_elementos
    na = conectividad_n(i, 1);
    nb = conectividad_n(i, 2);
    nodo_a = c_nod(na, :);
    nodo_b = c_nod(nb, :);
    
    % Longitud de la barra
    Lei = norm(nodo_b - nodo_a);
    L_barras(i) = Lei; 
    
    % Cosenos directores
    c(i) = (nodo_b(1) - nodo_a(1)) / Lei;
    s(i) = (nodo_b(2) - nodo_a(2)) / Lei;
    
    % Matriz de rotación nodal y elemental
    Re{i} = [c(i), s(i); 
          -s(i), c(i)];
    Te = [Re{i}, zeros(2); 
          zeros(2), Re{i}];
    
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

%% 6. Resolución del Sistema (desplazamientos)
gdl_libres = [7, 8, 9, 10];  % Grados de libertad libres (nodos 4 y 5)

K_LL = K(gdl_libres, gdl_libres);

% -------------------- Verificación de si es un mecanismo y restricciones en ejes locales del nodo ------------------------
if rank(K) < n_grados_libertad
    disp('La estructura es un mecanismo. No se puede resolver el sistema.');
    
    gdl_a_restringir = setdiff(gdl_libres, gdl_P); % Grados de libertad libres que no tiene la fuerza aplicada
    nodo_a_restringir = unique(ceil(gdl_a_restringir/2)); % Nodos a restringir
    barras_a_restringir = find(any(conectividad_n == nodo_a_restringir, 2)); % Barras conectadas a esos nodos

    % Display ordenado y tabulado
    fprintf('\n-------- ANÁLISIS DE RESTRICCIONES --------\n');
    fprintf('  GDL a restringir    : %s\n', num2str(gdl_a_restringir));
    fprintf('  Nodos a restringir  : %s\n', num2str(nodo_a_restringir));
    
    % Ponemos (:) para asegurar que se imprima en horizontal por si es vector columna
    fprintf('  Barras afectadas    : %s\n', num2str(barras_a_restringir(:)')); 
    fprintf('-------------------------------------------\n\n');

    theta = atan2(c(barras_a_restringir(2)), -s(barras_a_restringir(2))); % Ángulo en radianes entre los ejes globales y locales de la barra

    R = [cos(theta), sin(theta); -sin(theta), cos(theta)]; % Matriz de rotación para el nodo a restringir
    T = eye(n_grados_libertad); % Inicializamos la matriz de transformación global
    T(2*nodo_a_restringir-1:2*nodo_a_restringir, 2*nodo_a_restringir-1:2*nodo_a_restringir) = R; % T(9:10,9:10) Insertamos la matriz de rotación en la posición correspondiente

    Kt = T' * K * T; % Matriz de rigidez transformada
    Ft = T' * F; % Vector de fuerzas transformado

    gdl_lib_res = setdiff(gdl_libres, 2*nodo_a_restringir); % Grados de libertad libres después de la restricción (eliminamos el GDL restringido)

    u_l = zeros(n_grados_libertad, 1); % Inicializamos el vector de desplazamientos locales
    u_l(gdl_lib_res) = Kt(gdl_lib_res, gdl_lib_res) \ Ft(gdl_lib_res); % Calculamos los desplazamientos locales en los GDL libres restantes
    u_global = T * u_l; % Transformamos los desplazamientos locales a globales

else 
    disp('La estructura es estable. Se puede resolver el sistema.');
    % Cálculo de desplazamientos, vector de desplazamientos locales (F = K*u -> u = K\F)
    u_l = K_LL \ F(gdl_libres);

    % Reconstrucción del vector de desplazamientos global completo
    u_global = zeros(n_grados_libertad, 1);
    u_global(gdl_libres) = u_l;
end

% Cálculo de las nuevas coordenadas (deformadas)
c_nod_def = c_nod;
for i = 1:size(c_nod, 1)
    c_nod_def(i, :) = c_nod(i, :) + [u_global(2*i - 1), u_global(2*i)];
end

%% 7. Cálculo de tensiones en cada barra

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

%% 8. Representación Gráfica
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

title('Deformación de la Estructura');
xlabel('Posición X [m]');
ylabel('Posición Y [m]');
legend('Original', 'Deformada', 'Location', 'best');
hold off;