%% Problema 0 de cálculo matricial enb MEF
clear; close all; clc

%% Datos de entrada
L = 1; % longitud característica de la estructura
E = 1; % modulo elasticidad del material
A = 1; % area de las barras
I = 1; % momento de inercia de la viga continua
P = 1; % carga externa
gamma = pi/4; % ángulo de la carga sobre la estructura

% datos realistas
% L = 1; % longitud característica de la estructura
% E = 70e9; % modulo elasticidad del material
% A = 1e-4; % area de las barras
% I = 1e-8; % momento de inercia de la viga continua
% P = 1e4; % carga externa
 
%% Geometría y malla
n_1 = [0 L]; %coordenadas del nodo 1
n_2 = [0 0];
n_3 = [0 -1.5*L];
n_4 = [L 0];
n_5 = [L/2 L/2]; % NODO NUEVO

% operaciones con la malla de nodos:
c_nod = [n_1;n_2;n_3;n_4;n_5]; % matriz de coordenadas nodales
conectividad_n = [1 5;2 4;3 4;5 4]; % matriz de conectividad por nodos: nodos de cada elemento
nnod = size(c_nod,1); % número de nodos de la estructura
gdlnod = 2; % grados de libertad por nodo: en un elemento barra son 2
nGDL = nnod*gdlnod; % número de grados de libertad totales de la estructura: 2 gdls por nodo (elemento barra en el plano 2D)
nele = size(conectividad_n,1); % número de elementos de la estructura: el número de filas de conectividad_n

% graficamos estructura sin deformar
figure(1)
for e=1:nele
    nodos_e = conectividad_n(e,:); % nodos de los elementos
    n1 = c_nod(nodos_e(1),:); 
    n2 = c_nod(nodos_e(2),:);
    c_X = [n1(1) n2(1)]; 
    c_Y = [n1(2) n2(2)];
    plot(c_X,c_Y,'k-o')
    hold on
    grid on
    axis equal
end
title('estructura sin deformar y deformada')

%% Construcción de los elementos
% hoy ensamblamos en automático
% a partir de ahora: TODO en BUCLE AUTOMATICO
% Rigidez
K = zeros(nGDL); % inicio la matriz de rigidez vacía

% Fuerzas nodales equivalentes
Fneq = zeros(nGDL,1); % inicio vector de fuerzas nodales equivalentes de la estructura
conectividad_g = zeros(nele,4); % matriz de conectividad por grados de libertad: filas-elementos, columnas-gdls

for e=1:nele
    nodos = conectividad_n(e,:); % nodos de los elementos
    n1 = c_nod(nodos(1),:); % coordenadas del nodo local 1
    n2 = c_nod(nodos(2),:); % coordenadas del nodo local 2
    
    % grados de libertad de los elementos: me invento una regla nemotécnia
    % que funcione a partir de los nodos usados (o la busco por internet)
    gdle = [nodos(1)*2-1 nodos(1)*2 nodos(2)*2-1 nodos(2)*2];
    conectividad_g(e,:) = gdle; % guardos los grados de libertad
    
    [K_e,Te,Le,C,S] = rigidez_ele_axil(n1,n2,A,E); % rigidez del elemento barra
    K(gdle,gdle) = K(gdle,gdle) + K_e; % ensamblaje de la rigidez
    
    fneq_e = zeros(4,1); % fuerzas nodales equivalentes: no hay en este ejercicio
    Fneq(gdle) = Fneq(gdle) + fneq_e; % ensamblaje de las fuerzas nodales equivalentes
end

%% Fuerzas de acción en la estructura
Fn = zeros(nGDL,1); % vector de fuerzas nodales directas
gdl_P = [7 8]; % grados de libertad de la carga P
Fn(gdl_P) = P*[cos(gamma);-sin(gamma)]; % actualización de P

% vector de de fuerzas de acción
F = Fn + Fneq; % fuerzas de accion totales en los nodos

%% Condiciones de contorno
gdl_T = [1:nGDL]; % grados de libertad de la estructura
gdl_R = [1 2 3 4 5 6]; % grados restringidos
% gdl_L = [7 8 9 10]; % grados libres
gdl_L = setdiff(gdl_T,gdl_R);

% Condensación estática de matrices
U = zeros(nGDL,1); % vector de desplazamientos vacío
uR = U(gdl_R); % vector de desplazamientos restringidos
KLL = K(gdl_L,gdl_L); % matriz reducida por las C.C.: grados de deformación
KRR = K(gdl_R,gdl_R); % matriz reducida por las C.C.: grados restringidos
KLR = K(gdl_L,gdl_R); % matriz reducida por las C.C.: acoplamiento L-R
KRL = KLR'; % matriz reducida por las C.C.: acoplamiento L-R
fL = F(gdl_L); % fuerzas en grados libres: fuerzas de acción (dato)

%% Cálculo de incógnitas
uL = KLL^-1*(fL - KLR*uR); % desplazamiento en grados libres: incógnita
fR = KRL*uL + KRR*uR; % fuerzas externas en grados libres: incógnita
FR = fR - Fneq(gdl_R); % fuerzas de reacción en apoyos
F(gdl_R) = fR; % vector completo de fuerzas externas

% operaciones accesorias
U(gdl_L) = uL; % actualizo uL en U total
Ux = U([1:2:nGDL]); % desplazamientos nodales en eje x
Uy = U([2:2:nGDL]); % desplazamientos nodales en eje y
U_XY = [Ux Uy]; % matriz de desplazamientos nodales

% % U_ele = zeros(nele,4)
% for e=1:nele
%     gdle = conectividad_g(e,:);
%     uele = U(gdle);
%     U_ele(e) = uele; % desplazamientos nodales por elementos
% end

% Comprobación del equilibrio
Fint = K*U; % fuerzas internas en los nodos
Fext = F; % fuerzas externas en los nodos
error_equlibrio = Fint - Fext; % error en el equilibrio de fuerzas

if any(abs(error_equlibrio) > 1e-6)
    disp('error en el equilibrio final')
else
    disp('equilibrio correcto')
end

%% Representación gráfica de la estructura deformada
escala = L/(max(abs(U)))*0.3; % factor de escala para la representación gráfica
c_nod_def = c_nod + U_XY*escala; % coordenadas nodales deformadas
figure(1)
for e=1:nele
    nodos_e = conectividad_n(e,:); % nodos de los elementos
    n1 = c_nod_def(nodos_e(1),:); 
    n2 = c_nod_def(nodos_e(2),:);
    c_X = [n1(1) n2(1)]; 
    c_Y = [n1(2) n2(2)];
    plot(c_X,c_Y,'r--*')
    hold on
    grid on
    axis equal
end

%% Representación gráfica de la matriz de rigidez K
figure(2);
verK = zeros(nGDL);
for s=1:size(K,1)
    for t=1:size(K,1)
        if abs(K(s,t))>1e-12
            verK(s,t)=1;
        else 
            verK(s,t)=0; % Corregido: antes era una expresión lógica sin asignación
        end
    end
end
imagesc(verK); % Crea un plot de colores con los valores de K
colormap(flipud(gray)); % cambia el mapa de colores a gris
% los mayores negro, y los menores blanco
title('Ocupación de K','Fontsize',10)
xlabel('gdln','Fontsize',10) % Etiqueta el eje horizontal
ylabel('gdln','Fontsize',10) % Etiqueta el eje vertical

%% Cálculo de autovalores del problema: modos de deformación de la estructura con condiciones de contorno aplicadas
% CALCULO con KLL ahora, como hace NASTRAN. no uso la K completa: ya se lo que sale del ejercicio anterior
[autovecs,autovals] = eig(KLL); % autovectores y autovalores de la matriz KLL
[autovalores,ordeneig] = sort(diag(real(autovals)),'descend'); % ordenamiento de autovalores de mayor a menor
% help sort, help diag, help real
eigdim = length(autovalores); % el número de autovalores es el de grados de libertad totales ngdl
autovectores = zeros(length(gdl_L),eigdim);

% energía potencial de deformación
for m = 1:eigdim
    autovectores(:,m) = autovecs(:,ordeneig(m)); % ordenamiento de autovectores con autovalores de mayor a menor
    % energía potencial interna de deformación del modo "m"
    PIint_e(m) = 1/2*autovectores(:,m)'*KLL*autovectores(:,m); % ESTE era el FALLO mío en CLASE: no usar KLL, sino K
end

%% Representa gráficamente los modos de deformación y de sólido rígido
% Autovalores y Modos de deformación
figure(3)
mm = ceil(sqrt(length(autovalores))); % entero superior de sqrt()
figure(3); set(figure(3),'WindowState','maximize');

for n = 1:eigdim
    autovalor = autovalores(n);
    U_eig = zeros(nGDL,1); % declaro el tamaño del vector de desplazamientos nodales vacío
    U_eig(gdl_L) = real(autovectores(:,n)); % añado los desplazamientos de gdl_L: los de KLL
    escala_eig = 0.1*L/(max(abs(U_eig))); % escala automática de representación
    subplot(mm,mm,n);
    
    for e=1:nele
        gdle = conectividad_g(e,:);
        node = conectividad_n(e,:);
        cne = [c_nod(node(1),:),c_nod(node(2),:)];
        cx_n = cne([1 3]); cy_n = cne([2 4]);
        plot(cx_n,cy_n,'ko-','LineWidth',1)
        hold on
        ue = real(U_eig(gdle)); 
        uex = ue([1 3])'; 
        uey = ue([2 4])';
        cx_eig = cx_n+escala_eig*uex; 
        cy_eig = cy_n+escala_eig*uey;
        cplot = plot(cx_eig,cy_eig,'ro--','LineWidth',1);
        % axis equal
    end
    
    xlabel('coordenadas X','Fontsize',12)
    ylabel('coordenadas Y','Fontsize',12)
    % zticklabels('z = 0')
    grid on
    axis equal
    format shortE
    title(['Modo ' int2str(n) ' de{\it K},' ' autovalor ' num2str(autovalor, '%6.3e') ','...
    sprintf('\n') 'energía de deformación = ' num2str(PIint_e(n), '%6.3e')],'Fontsize',12)
end
% saveas(gcf,'modos_deformacion.png')

%% Cálculo de tensiones: en MEF ya...

%% Subrutinas: "functions" de Matlab
function [K_e,Te,Le, C,S] = rigidez_ele_axil(n1,n2,A,E)
    % n1, n2: coordenadas globales de los nodos locales 1 y 2 del elemento de 2 nodos
    % E: módulo de Young
    % A: área de sección recta del elemento - constante
    
    % distancias relativas entre nodos locales 1 y 2 en el plano 2D
    deltax_e = n2(1)-n1(1); 
    deltay_e = n2(2)-n1(2);
    Le = sqrt(deltax_e^2+deltay_e^2);
    
    % ángulo entre coordenadas locales y globales en el plano 2D
    S = deltay_e/Le; % seno de alfa
    C = deltax_e/Le; % coseno de alfa
    
    % vectores de la base local xy en coordenadas globales
    n_x = [1;0;0]; 
    n_y = [0;1;0]; 
    n_z = [0;0;1]; % base de referencia global
    
    % vectores de la base local x'y' en coordenadas globales
    nx = [C;S;0]; 
    ny = [-S;C;0]; 
    nz = [0;0;1];
    
    % matriz de cambio de base de un vector 2D
    % Re = [nx'*n_x nx'*n_y nx'*n_z;ny'*n_x ny'*n_y ny'*n_z;nz'*n_x nz'*n_y nz'*n_z]; %Teoría
    Re = [C S;-S C]; %alternativa acortada
    
    %desacoplamiento de los nodos del elemento entre si
    Zeros = zeros(2);
    
    % matriz de transformación de un elemento de 2 nodos en el plano 2D
    Te = [Re Zeros;Zeros Re];
    Ke = A*E/Le*[1 0 -1 0;0 0 0 0;-1 0 1 0;0 0 0 0]; % matriz de rigidez del elemento en coordenadas locales
    K_e = Te'*Ke*Te; % matriz de rigidez del elemento en coordenados globales
end