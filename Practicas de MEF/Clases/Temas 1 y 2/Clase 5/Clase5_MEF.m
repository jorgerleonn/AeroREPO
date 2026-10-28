%% Vigas continuas con rótula en MEF
clear; close all; fn = 1;
% Valores realistas:
L = 1000; % mm
a = L/20; A = a^2; I = a^4/12; % sección recta de viga cuadrada
E = 100e3; % N mm^-2
P1 = -10e3; % N
gdl_P1 = 5; % grado de libertad donde se halla P1
q0y = -10; % fuerza local en el interior de los elementos N mm^-1

%% Geometría y Malla: nodos y elementos
% Sistema de referencia: direcciones de la base global en la matriz DIR
dirX = [1 0 0]'; dirY = [0 1 0]'; dirZ = [0 0 1]'; DIR = [dirX,dirY,dirZ];
% Nodos escritos a mano
n_1 = [0 0];
n_2 = [L/2 0];
n_3 = [L 0];
n_4 = [L 0];
n_5 = [1.5*L 0];
n_6 = [2*L 0];
c_nod = [n_1;n_2;n_3;n_4;n_5;n_6]; % matriz de coordenadas nodales
nnod = size(c_nod,1); % número de nodos
% Elementos: matriz de conectividad por nodos
% elementos en filas y nodos en columnas 1 y 2
conectividad_n = [1 2;2 3;4 5;5 6];
nele = size(conectividad_n,1); % número de elementos de la estructura
% Grados de libertad
gdln = 3; % grados de libertad por nodo
nGDL = gdln*nnod; % grados de libertad totales

%% Cálculo dimensional-espacial de los elementos
L_e = zeros(1,nele); % almacen de dimensiones
T_e = zeros(6,6,nele); % almacen de matrices de cambio de base
conectividad_g = zeros(nele,6); % conectividad_g automática
for e=1:nele
    nodos = conectividad_n(e,:);
    n1 = c_nod(nodos(1),:); n2 = c_nod(nodos(2),:);
    [Le,Te] = Geom_viga2D(n1,n2,DIR);
    L_e(e) = Le;
    T_e(:,:,e) = Te;
    % conectividad_g automática
    conectividad_g(e,:) = [nodos(1)*gdln-2 nodos(1)*gdln-1 nodos(1)*gdln ...
                           nodos(2)*gdln-2 nodos(2)*gdln-1 nodos(2)*3];
end

%% Fuerzas locales en los elementos (Condiciones de contorno en fuerzas: naturales o de Neumann)
% vector de fuerzas nodales equivalentes de los elementos en coordenadas locales
Fneq_ele = zeros(6,nele); % almacen de fuerzas nodales equivalentes por elementos-columnas
for e=1:nele % método en bucle
    Le = L_e(e);
    if e==1 % elemento 1
        fneqe = q0y*1/2*Le*[0;3/20;Le/30;0;7/20;-Le/20];
    elseif e==2 % elemento 2
        fneqe = q0y*1/2*Le*[0;1/2;Le/12;0;1/2;-Le/12] + ...
                q0y*1/2*Le*[0;3/20;Le/30;0;7/20;-Le/20];
    elseif e==3 % elemento 3
        fneqe = q0y*Le*[0;1/2;Le/12;0;1/2;-Le/12];
    elseif e==4 % elemento 4
        fneqe = q0y*Le*[0;1/2;Le/12;0;1/2;-Le/12];
    else
        fneqe = [0;0;0;0;0;0];
    end
    Fneq_ele(:,e) = fneqe; % almacenamiento
end

%% ensamblaje en la matriz de rigidez global
K = zeros(nGDL); % matriz de rigidez global de la estructura
K_ele = zeros(6,6,nele);
Fneq = zeros(nGDL,1);
for e=1:nele
    gdl_e = conectividad_g(e,:);
    Le = L_e(e);
    Ke = Ke_viga2D(E,Le,A,I); % rigidez local
    Te = T_e(:,:,e);
    K_e = Te'*Ke*Te; % rigidez global
    K(gdl_e,gdl_e) = K(gdl_e,gdl_e) + K_e; % ensamblaje de rigidez
    K_ele(:,:,e) = K_e;
    % fuerzas nodales equivalentes
    fneqe = Fneq_ele(:,e); % fuerzas neq locales
    fneq_e = Te'*fneqe; % fuerzas neq globales
    Fneq(gdl_e) = Fneq(gdl_e) + fneq_e; % ensamblaje de fuerzas
end
rangoK = rank(K);
deficit_rangoK = size(K,1)-rangoK;

%% Fuerzas globales de acción en los nodos: nodales directas
Fn = zeros(nGDL,1); Fn(gdl_P1) = P1;

%% Fuerzas externas nodales de acción
FA = Fn + Fneq; % sin condiciones de contorno

%% Condiciones de contorno en desplazamientos (esenciales o de Dirichlet)
U = zeros(nGDL,1); % vector inicial de desplazamientos
F = zeros(nGDL,1);
gdlT = [1:nGDL]; % grados de libertad totales
gdlR = [1 2 3 nGDL-1]; % grados de libertad restringidos: SPC
uR = [0;0;0;0]; % vector de desplazamientos impuestos: SPC
% U(gdlR) = uR; % sustitución de desplazamientos impuestos con rot incluida
gdlL = setdiff(gdlT,gdlR);

%% Solución con matriz de coeficientes de restricción
% SPC
r_SPC = uR;
p_SPC = length(uR); % numéro de restricciones
% r = zeros(p_SPC,1); % restricciones totales de la estructura
C_SPC = zeros(p_SPC,nGDL); % matriz de coeficientes de restricción o de constraints
C_SPC(1,1) = 1; % u1 = 0
C_SPC(2,2) = 1; % u2 = 0
C_SPC(3,3) = 1; % u3 = 0
C_SPC(4,17) = 1; % u17 = 0

% MPC
p_MPC = 2; %
r_MPC = zeros(p_MPC,1); % restricciones totales de la estructura
C_MPC = zeros(p_MPC,nGDL); % matriz de coeficientes de restricción o de constraints
C_MPC(1,7) = 1; C_MPC(1,10) = -1; % u7 - u10 = 0
C_MPC(2,8) = 1; C_MPC(2,11) = -1; % u8 - u11 = 0

% Total
p_Total = p_SPC + p_MPC;
r = [r_SPC;r_MPC];
C = [C_SPC;C_MPC]; % matriz total de coeficientes de restricción

%% Multiplicadores de Lagrange - Verificación
K_ML = [K C';C zeros(p_Total)]; % matriz aumentada de los ML
F_ML = [FA;r]; % vector de "fuerzas" ampliado por los ML
U_ML = K_ML\F_ML; % vector de incógnitas ML

% incógnitas
U_ml = U_ML(gdlT); % vector de desplazamientos incógnita ML
lambda = U_ML([length(gdlT)+1:end]);
FR_ml = -lambda; % fuerzas de reacción condensadas
FR_ML = C'*FR_ml; % fuerzas de reacción totales sin condensar

FR = FR_ml([1:4]); % fuerzas de reacción en apoyos

F_rotula_I = FR_ML([7 8]); % fuerzas internas en la rótula
F_rotula_D = FR_ML([10 11]); % fuerzas internas en la rótula

%% Verificación ML
Fint = K*U_ml;
Fext = FA + FR_ML;
errorF_ml = Fint - Fext; % verificación de equilibrio interno
if all(abs(errorF_ml) < 1e-6)
    disp('Funciona OK')
else
    disp('error: no funciona')
end

%% Método de la penalización
kpen = 1e4*max(max(K)); % parámetro escalar de penalizacion
k_pen = kpen*eye(p_Total); % matriz de coeficientes de penalizacion
K_Pen = C'*k_pen*C; % matriz completa de penalización
F_Pen = C'*k_pen*r; % fuerzas de penalización

% incógnitas
U_Pen = (K + K_Pen)\(FA + F_Pen); % vector de desplazamientos penalizados
Fint = K*U_Pen;
FR_Pen = Fint - FA;

%% Verificación Pen
U = U_Pen;
Fint = K*U;
Fext = FA + FR_Pen;
errorF_pen = Fint - Fext; % verificación de equilibrio interno
if all(abs(errorF_pen) < 1e-6)
    disp('Funciona OK')
else
    disp('error: no funciona')
end

%% Diferencias entre métodos
dif_U = U_ml - U_Pen;
dif_FR = FR_ML - FR_Pen;

%% Deformada
U = U_ml;
Ux = U([1:3:nGDL]);
Uy = U([2:3:nGDL]);
Uz = U([3:3:nGDL]);
c_nod_def = c_nod + [Ux,Uy];

%% Fuerzas y desplazamientos por elementos
U_e = zeros(6,nele);
F_e = zeros(6,nele);
Fint_e = zeros(6,nele);
for e=1:nele
    gdl_e = conectividad_g(e,:);
    u_e = U(gdl_e); % desplazamientos en globales del ele1
    U_e(:,e) = u_e; % desplazamientos de los elementos
    Fint_e(:,e) = (K_ele(:,:,e)*u_e); % fuerzas nodales internas
end

%% Representación gráfica
escala = L/max(abs([Ux;Uy]))*0.1; % escala de magnificación 10%
cnoddef_esc = c_nod + [Ux,Uy]*escala; % coordenadas de los nodos deformados escalados

figure(fn); fn = fn + 1;
plot(c_nod(:,1),c_nod(:,2),'b-o')
hold on
grid on
axis equal
title('Estructura original y deformada elástica','Fontsize',12)
plot(cnoddef_esc(:,1),cnoddef_esc(:,2),'r--*')

figure(fn); fn = fn + 1;
verK = zeros(nGDL);
for s=1:size(K,1)
    for t=1:size(K,1)
        if abs(K(s,t))>1e-12
            verK(s,t)=1;
        end
    end
end
imagesc(verK); % Crea un plot de colores con los valores de K
colormap(flipud(gray)); % cambia el mapa de colores a gris
title('Ocupación de K','Fontsize',10)
xlabel('gdln','Fontsize',10) % Etiqueta el eje horizontal
ylabel('gdln','Fontsize',10) % Etiqueta el eje vertical

%% Subrutinas: "functions"
function [Le,Te] = Geom_viga2D(n1,n2,DIR)
% n1,n2: coordenadas globales de los nodos locales 1 y 2
% base global de referencia
n_X = DIR(:,1); n_Y = DIR(:,2); n_Z = DIR(:,3);
% diferencias de longitud en coordenadas globales
delta = n2-n1; delta_x = delta(1); delta_y = delta(2); % distancia entre nodos
Le = sqrt(delta_x^2+delta_y^2); % Le longitud global del elemento
C = delta_x/Le; % coseno del ángulo en el plano
S = delta_y/Le; % seno del ángulo en el plano
% vectores de la base local x'y' en coordenadas globales
n_x = [C;S;0]; n_y = [-S;C;0]; n_z = [0;0;1];
% Re matriz de cambio de base global-local
CXx = n_x'*n_X; CYx = n_x'*n_Y; CZx = n_x'*n_Z;
CXy = n_y'*n_X; CYy = n_y'*n_Y; CZy = n_y'*n_Z;
CXz = n_z'*n_X; CYz = n_z'*n_Y; CZz = n_z'*n_Z;
Re = [CXx CYx CZx;CXy CYy CZy;CXz CYz CZz]; % Te matriz de transformación vector
R0 = zeros(3);
Te = [Re R0;R0 Re]; % Te matriz de transformación global-local de elemento de 2 nodos
end

function [Ke] = Ke_viga2D(E,Le,A,I)
% rigideces escalares
% E: módulo de Young
% A: área de sección recta del elemento - constante
% I: momento de inercia de la sección recta del elemento plano - constante
EA = E*A; EI = E*I;
k11 = [EA/Le 0 0; 0 12*EI/Le^3 6*EI/Le^2; 0 6*EI/Le^2 4*EI/Le];
k12 = [-EA/Le 0 0; 0 -12*EI/Le^3 6*EI/Le^2; 0 -6*EI/Le^2 2*EI/Le];
k22 = [EA/Le 0 0; 0 12*EI/Le^3 -6*EI/Le^2; 0 -6*EI/Le^2 4*EI/Le];
k21 = k12';
Ke = [k11 k12; k21 k22];
end