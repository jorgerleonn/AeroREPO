%% Clase 4 - ejercicio de Vigas continuas en MEF de las diapositivas
clear; close all;

%% 1. Parámetros de la estructura y del material
L = 1000; % mm
a = L/20; A = a^2; I = a^4/12; % sección recta de viga cuadrada
E = 100e3; % N mm^-2

P1 = -0e3; % N
gdl_P1 = 8; % grado de libertad donde se halla P1
P2 = -100e3; % N
gdl_P2 = 14; % grado de libertad donde se halla P2
q0y = -10; % fuerza local en el interior de los elementos N mm^-1

%% 2. Geometría y Malla: nodos y elementos
% Sistema de referencia: direcciones de la base global en la matriz DIR
dirX = [1 0 0]'; dirY = [0 1 0]'; dirZ = [0 0 1]'; DIR = [dirX, dirY, dirZ];

% Nodos y grados de libertad
n1 = [0 0];
n2 = [L/2 L/2];
n3 = [L L];
n4 = [1.5*L L];
n5 = [2*L L];
n6 = [2*L 0.5*L];
n7 = [2*L 0];
c_nod = [n1; n2; n3; n4; n5; n6; n7]; % matriz de coordenadas nodales
nnod = size(c_nod,1); % número de nodos de la estructura
% Elementos y conectividad por nodos
conectividad_n = [1 2; 2 3; 3 4; 4 5; 5 6; 6 7]; % matriz de conectividad por nodos (elementos en filas y nodos en columnas 1 y 2)
n_elementos = size(conectividad_n,1); % número de elementos de la estructura

gdlxnod = 3; % Grados de libertad por nodo (2D)
n_grados_libertad = gdlxnod * nnod; % Número de grados de libertad

%% 3. Cálculo dimensional-espacial de los elementos
Le = zeros(n_elementos, 1); % Longitud de los elementos
Te = zeros(6, 6, n_elementos); % Matriz de transformación de los elementos
conectividad_g = zeros(n_elementos, gdlxnod*2); % Matriz de conectividad por grados de libertad

for i = 1:n_elementos
    nodos = conectividad_n(i,:);
    na = c_nod(nodos(1),:); 
    nb = c_nod(nodos(2),:);
    [Le_i, Te_i] = Geom_viga2D(na, nb, DIR); % Longitud y matriz de transformación del elemento
    Le(i) = Le_i; % Guardar longitud del elemento
    Te(:,:,i) = Te_i; % Guardar matriz de transformación del elemento

    gdl_i = [gdlxnod*nodos(1)-2, gdlxnod*nodos(1)-1, gdlxnod*nodos(1), gdlxnod*nodos(2)-2, gdlxnod*nodos(2)-1, gdlxnod*nodos(2)]; % Grados de libertad del elemento
    conectividad_g(i, :) = gdl_i; % Guardar conectividad por grados de libertad
end

%% 5. Ensamblaje de la matriz de rigidez global y vector de fuerzas nodales equivalentes globales
K = zeros(n_grados_libertad, n_grados_libertad); % Matriz de rigidez global
K_ele = zeros(6,6, n_elementos); % Matriz de rigidez de los elementos
Fneq_ele = zeros(6, n_elementos);
Fneq = zeros(n_grados_libertad, 1); % Vector de fuerzas nodales equivalentes globales

for i = 1:n_elementos
    % Cálculo del vector de fuerzas nodales equivalentes globales
    if i == 1
        fneq = q0y*Le(i)*[0; 1/2; Le(i)/12; 0; 1/2; -Le(i)/12]; % fneqe = [u1; v1; θ1; u2; v2; θ2]
    elseif i == 2
        fneq = q0y*Le(i)*[0; 1/2; Le(i)/12; 0; 1/2; -Le(i)/12]; % fneqe = [u1; v1; θ1; u2; v2; θ2]
    else
        fneq = [0; 0; 0; 0; 0; 0];
    end

    % Ensamblaje de la matriz de rigidez global y vector de fuerzas nodales equivalentes globales
    gdl_i = conectividad_g(i, :); % Grados de libertad del elemento i
    Te_i = Te(:,:,i); % Matriz de transformación del elemento i
    Ke_i = Ke_viga2D(E, Le(i), A, I); % Matriz de rigidez local del elemento i
    K_e = Te_i'*Ke_i*Te_i; % Matriz de rigidez global del elemento i
    K(gdl_i, gdl_i) = K(gdl_i, gdl_i) + K_e; % Ensamblaje de la matriz de rigidez global
    K_ele(:,:,i) = K_e; % Guardar matriz de rigidez del elemento i

    Fneq(gdl_i) = Fneq(gdl_i) + Te_i'*fneq;
end

Fn = zeros(n_grados_libertad, 1); % Vector de fuerzas nodales equivalentes globales
Fn(gdl_P1) = P1; % Aplicar fuerza P1 en el grado de libertad correspondiente
Fn(gdl_P2) = P2; % Aplicar fuerza P2 en el grado de libertad correspondiente

F_total = Fn + Fneq; % Vector de fuerzas nodales equivalentes globales total; 

%% 6. Condiciones de contorno y Resolución del Sistema (desplazamientos)
U = zeros(n_grados_libertad, 1); % Vector de desplazamientos nodales globales
gdlT = [1:n_grados_libertad]; % Grados de libertad totales
gdlR = [1, 2, 3, 19, 20, 21]; % Grados de libertad restringidos
gdl_libres = setdiff(gdlT, gdlR); % Grados de libertad libres

K_LL = K(gdl_libres, gdl_libres); % Submatriz de rigidez global correspondiente a los grados de libertad libres
K_LR = K(gdl_libres, gdlR); % Submatriz de rigidez global correspondiente a los grados de libertad libres y restringidos
K_RR = K(gdlR, gdlR); % Submatriz de rigidez global correspondiente a los grados de libertad restringidos
K_RL = K_LR'; % Submatriz de rigidez global correspondiente a los grados de libertad restringidos y libres

u_R = zeros(length(gdlR), 1); % Desplazamientos restringidos
u_l = K_LL \ F_total(gdl_libres); % Resolver el sistema de ecuaciones para los desplazamientos libres
U(gdlR) = u_R; % Asignar desplazamientos restringidos al vector de desplazamientos nodales globales
U(gdl_libres) = u_l; % Asignar desplazamientos libres al vector de desplazamientos nodales globales

Ux = U(1:3:end); % Desplazamientos en x
Uy = U(2:3:end); % Desplazamientos en y
Uz = U(3:3:end); % Desplazamientos en z

% Coordenadas de los nodos deformados
c_nod_def = zeros(nnod, 2); % Inicializar matriz de coordenadas de los nodos deformados
for i = 1:nnod
    c_nod_def(i, 1) = c_nod(i, 1) + Ux(i); % Coordenada x del nodo deformado
    c_nod_def(i, 2) = c_nod(i, 2) + Uy(i); % Coordenada y del nodo deformado
end

% Fuerzas totales
F_reaccion = F_total(gdlR) - Fneq(gdlR); % Fuerzas de reacción en los grados de libertad restringidos

% Verificación de la solución
F_int = K * U; % Fuerzas internas
F_ext = F_total; % Fuerzas externas
error_F = F_int - F_ext; % Error en las fuerzas
if any(abs(error_F) > 1e-6)
    warning('Error en las fuerzas: %e', max(abs(error_F)));
else
    disp('Solución verificada: fuerzas internas y externas coinciden.');
end

%% 7. Fuerzas y desplazamientos por elementos
U_e =zeros(6, n_elementos); % Desplazamientos por elementos
F_e = zeros(6, n_elementos); % Fuerzas por elementos
Fint_e = zeros(6, n_elementos); % Fuerzas internas por elementos

for i = 1:n_elementos
    gdl_i = conectividad_g(i, :); % Grados de libertad del elemento i
    u_e = U(gdl_i); % Desplazamientos del elemento i
    U_e(:, i) = u_e; % Guardar desplazamientos del elemento i
    Fint_e(:, i) = K_ele(:, :, i) * u_e; % Calcular fuerzas internas del elemento i
end

%% 8. Autovalores y modos de deformación
[autovecs, autovals] = eig(K);
[autovalores, orden_a] = sort(diag(real(autovals)), 'descend'); % Ordenamos los autovalores

n_modos = length(orden_a); % Número de modos
autovectores = zeros(n_grados_libertad, n_modos); % Inicializamos la matriz de autovectores
PIint_e = zeros(1, n_modos); % Inicializamos el vector de energía de deformación interna por modo

for i = 1:n_modos
    autovectores(:, i) = autovecs(:, orden_a(i)); % Reordenamos los autovectores según los autovalores
    PIint_e(i) = 0.5 * autovectores(:, i)' * K * autovectores(:, i); % Energía de deformación interna por modo
end

%% Representaciones gráficas
% Representación gráfica de la estructura sin deformar y deformada
figure(1);
for i = 1:n_elementos
    nodos = conectividad_n(i,:);
    na = c_nod(nodos(1),:); 
    nb = c_nod(nodos(2),:);
    na_def = c_nod_def(nodos(1),:);
    nb_def = c_nod_def(nodos(2),:);
    plot([na(1) nb(1)], [na(2) nb(2)], 'b-o')
     hold on;
    plot([na_def(1) nb_def(1)], [na_def(2) nb_def(2)], 'r--*')
    hold on;
end
grid on
axis equal
title('Estructura original y deformada elástica', 'Fontsize', 12)
hold off;

% Representación gráfica de la estructura deformada con escala de magnificación
figure(2);
escala = L/max(abs([Ux; Uy]))*0.1; % escala de magnificación 10%
cnod_def_esc = c_nod + [Ux, Uy]*escala; % coordenadas de los nodos deformados escalados
for i = 1:n_elementos
    nodos = conectividad_n(i,:);
    na = c_nod(nodos(1),:); 
    nb = c_nod(nodos(2),:);
    na_def = cnod_def_esc(nodos(1),:);
    nb_def = cnod_def_esc(nodos(2),:);
    plot([na(1) nb(1)], [na(2) nb(2)], 'b-o')
     hold on;
    plot([na_def(1) nb_def(1)], [na_def(2) nb_def(2)], 'r--*')
    hold on;
end
grid on
axis equal
title('Estructura original y deformada elástica', 'Fontsize', 12)
hold off;

% Representación de la ocupación de la matriz de rigidez global
figure(3);
verK = zeros(n_grados_libertad);
for i=1:size(K,1)
    for j=1:size(K,1)
        if abs(K(i,j)) > 1e-12
            verK(i,j) = 1;
        else
            verK(i,j) = 0;
        end
    end
end
imagesc(verK);
colormap(flipud(gray));
title('Ocupación de K','Fontsize',10);
xlabel('Grados de libertad','Fontsize',10);
ylabel('Grados de libertad','Fontsize',10);

% Representación de los modos de deformación y de sólido rígido
figure('Name','Modos de deformación de K','Color','w');
M = [1 2 10 15 18 19 20 21]; % Modos a representar (extremos del espectro + rígidos)
mm = ceil(sqrt(length(M))); % Número de subplots por fila y columna

for n = 1:length(M)
    m = M(n); % Modo a representar
    U_eig = real(autovectores(:,m)); % Desplazamientos del modo n
    escala_eig = 0.3*L / max(abs(U_eig)); % escala automática

    subplot(mm, mm, n); hold on; grid on; axis equal;
    for e = 1:n_elementos
        gdle    = conectividad_g(e,:);
        nodos_e = conectividad_n(e,:);

        cx_n = [c_nod(nodos_e(1),1), c_nod(nodos_e(2),1)];
        cy_n = [c_nod(nodos_e(1),2), c_nod(nodos_e(2),2)];
        plot(cx_n, cy_n, 'ko-','LineWidth',1);       % modo sin deformar

        ue  = real(U_eig(gdle));                     % 6 desplazamientos del elemento
        uex = ue([1 4])';                            % u en nodos 1 y 2
        uey = ue([2 5])';                            % v en nodos 1 y 2
        plot(cx_n + escala_eig*uex, cy_n + escala_eig*uey, 'ro--','LineWidth',1);
    end
    title(sprintf('Modo %d | autovalor = %.3e | energía = %.3e', ...
                  m, autovalores(m), PIint_e(m)));
    xlabel('X [m]'); ylabel('Y [m]');
end