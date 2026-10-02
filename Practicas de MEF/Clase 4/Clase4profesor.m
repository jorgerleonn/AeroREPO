%% Clase 4 ejercicio de Vigas continuas en MEF de las diapositivas
clear; close all; fn = 1;

% Valores realistas:
L = 1000; % mm
a = L/20; A = a^2; I = a^4/12; % sección recta de viga cuadrada
E = 100e3; % N mm^-2
P1 = -0e3; % N
gdl_P1 = 8; % grado de libertad donde se halla P1
P2 = -100e3; % N
gdl_P2 = 14; % grado de libertad donde se halla P2
q0y = -10; % fuerza local en el interior de los elementos N mm^-1

%% 1. Geometría y Malla: nodos y elementos
% Sistema de referencia: direcciones de la base global en la matriz DIR
dirX = [1 0 0]'; dirY = [0 1 0]'; dirZ = [0 0 1]'; DIR = [dirX, dirY, dirZ];

% Nodos escritos a mano
n_1 = [0 0];
n_2 = [L/2 L/2];
n_3 = [L L];
n_4 = [1.5*L L];
n_5 = [2*L L];
n_6 = [2*L 0.5*L];
n_7 = [2*L 0];
c_nod = [n_1; n_2; n_3; n_4; n_5; n_6; n_7]; % matriz de coordenadas nodales
nnod = size(c_nod,1); % número de nodos de la estructura

% Elementos: matriz de conectividad por nodos (elementos en filas y nodos en columnas 1 y 2)
conectividad_n = [1 2; 2 3; 3 4; 4 5; 5 6; 6 7];
nele = size(conectividad_n,1); % número de elementos de la estructura

% Representación gráfica de la estructura
% apaño rápido para elementos 1D
fn = 1; figure(fn); fn = fn + 1;
% plot(c_nod(:,1),c_nod(:,2),'b-o') % esto funciona solo si los elementos son continuos todos
% hold on
% grid on
% axis equal
% title('Estructura original y deformada elástica','Fontsize',12)
for e = 1:nele % esto funciona siempre
    conecta_nodos = conectividad_n(e,:);
    % elemento sin deformar
    n1 = c_nod(conecta_nodos(1),:); 
    n2 = c_nod(conecta_nodos(2),:);
    plot([n1(1) n2(1)], [n1(2) n2(2)], 'b-o')
    hold on
end
grid on
axis equal
title('Estructura original y deformada elástica', 'Fontsize', 12)

% Grados de libertad
gdln = 3; % grados de libertad por nodo
nGDL = gdln*nnod; %grados de libertad totales

% % Elementos: matriz de conectividad por grados de libertad
% gdle1 = [1 2 3 4 5 6]; % A manubrio
% gdle2 = [4 5 6 7 8 9];
% gdle3 = [7 8 9 10 11 12];
% gdle4 = [10 11 12 13 14 15];
% gdle5 = [13 14 15 16 17 18];
% gdle6 = [16 17 18 19 20 21];
% % elementos en filas y glds en columnas
% conectividad_g = [gdle1;gdle2;gdle3;gdle4;gdle5;gdle6]; % A manubrio

% cálculo dimensional-espacial de los elementos
L_e = zeros(1, nele); % almacen de dimensiones
T_e = zeros(6, 6, nele); % almacen de matrices de cambio de base
conectividad_g(nele, 6) = zeros; % conectividad_g automática

for e = 1:nele
    nodos = conectividad_n(e,:);
    n1 = c_nod(nodos(1),:); 
    n2 = c_nod(nodos(2),:);
    [Le, Te] = Geom_viga2D(n1, n2, DIR);
    L_e(e) = Le;
    T_e(:,:,e) = Te;
    % conectividad_g automática
    conectividad_g(e,:) = [nodos(1)*gdln-2, nodos(1)*gdln-1, nodos(1)*gdln, nodos(2)*gdln-2, nodos(2)*gdln-1, nodos(2)*3];
end

%% Fuerzas locales en los elementos (Condiciones de contorno en fuerzas: naturales o de Neumann)
% vector de fuerzas nodales equivalentes de los elementos en coordenadas locales
% a manubrio
% Le1 = L_e(1);
% fneqe1 = q0y*Le1*[0;1/2;Le1/12;0;1/2;-Le1/12];
% Le2 = L_e(2);
% fneqe2 = q0y*Le2*[0;1/2;Le2/12;0;1/2;-Le2/12];
% fneqe3 = [0;0;0;0;0;0];
% fneqe4 = [0;0;0;0;0;0];
% fneqe5 = [0;0;0;0;0;0];
% fneqe6 = [0;0;0;0;0;0];
% % almacen de fuerzas nodales equivalentes por elementos-columnas
% Fneq_ele = [fneqe1,fneqe2,fneqe3,fneqe4,fneqe5,fneqe6];

% automático
Fneq_ele = zeros(6, nele); % almacen de fuerzas nodales equivalentes por elementos-columnas
for e = 1:nele % método en bucle
    Le = L_e(e);
    if e == 1 % elemento 1
        fneqe = q0y*Le*[0; 1/2; Le/12; 0; 1/2; -Le/12];
    elseif e == 2 % elemento 2
        fneqe = q0y*Le*[0; 1/2; Le/12; 0; 1/2; -Le/12];
    else
        fneqe = [0; 0; 0; 0; 0; 0];
    end
    Fneq_ele(:,e) = fneqe; % almacenamiento
end

%% ensamblaje en la matriz de rigidez global
K = zeros(nGDL); % matriz de rigidez global de la estructura
K_ele = zeros(6, 6, nele);
Fneq = zeros(nGDL, 1);

for e = 1:nele
    gdl_e = conectividad_g(e,:);
    Le = L_e(e);
    Ke = Ke_viga2D(E, Le, A, I); % rigidez local
    Te = T_e(:,:,e);
    K_e = Te'*Ke*Te; % rigidez global
    K(gdl_e, gdl_e) = K(gdl_e, gdl_e) + K_e; %ensamblaje de rigidez
    K_ele(:,:,e) = K_e;
    
    % fuerzas nodales equivalentes
    fneqe = Fneq_ele(:,e); % fuerzas neq locales
    fneq_e = Te'*fneqe; % fuerzas neq globales
    Fneq(gdl_e) = Fneq(gdl_e) + fneq_e; %ensamblaje de fuerzas
end

%% Fuerzas globales de acción en los nodos: nodales directas
Fn = zeros(nGDL, 1);
Fn(gdl_P1) = P1;
Fn(gdl_P2) = P2;
% Fuerzas externas nodales de acción
FA = Fn + Fneq;

%% Autovalores y Modos de deformación
[eigenvecs, eigenvals] = eig(K); % autovectores y autovalores de la matriz K
[autovalores, orden_a] = sort(diag(real(eigenvals)), 'descend'); % ordenamiento de autovalores
autovectores = zeros(length(autovalores));
PI_int = zeros(1, length(orden_a)); % Inicializar vector para evitar warnings

for i = 1:length(orden_a)
    i_a = orden_a(i);
    autovectores(:,i) = eigenvecs(:,i_a); % ordenamiento de autovectores
    PI_int(i) = 1/2*autovectores(:,i)'*K*autovectores(:,i); % energia de deformación
end

PI_int([1 2 3]) % energia de los primeros 3 modos
PI_int([end-2 end-1 end]) % energía de los últimos 3 modos
M = [1 2 10 15 18 19 20 21]; % modos elegidos a representar
mm = ceil(sqrt(length(M)));

fign = figure(fn); fn = fn + 1;
set(fign, 'WindowState', 'maximize');

for n = 1:length(M)
    m = M(n);
    autovalor = autovalores(m);
    U_eigen = real(autovectores(:,m));
    
    % energía de deformación
    %PI_int(m) = 1/2*U_eigen'*K*U_eigen;
    % m = aindex(1);
    % U_eigen = autovecs(:,m);
    
    escala2 = 0.3*L/(max(abs(U_eigen)));
    subplot(2, 4, n);
    
    for e = 1:nele
        gdle = conectividad_g(e,:);
        nep = conectividad_n(e,:);
        cne = [c_nod(nep(1),:), c_nod(nep(2),:)];
        cx_n = cne([1 3]); 
        cy_n = cne([2 4]);
        plot(cx_n, cy_n, 'bo-', 'LineWidth', 2)
        hold on
        
        %ue = U(gdl_e(e,:));
        %ue = U_e(e,:)';
        ue = real(U_eigen(gdle));
        uex = ue([1 4])'; 
        uey = ue([2 5])';
        cplot = plot(cx_n+escala2*uex, cy_n+escala2*uey, 'ro--', 'LineWidth', 2);
        % axis equal
    end
    
    xlabel('coordenadas X', 'Fontsize', 12)
    ylabel('coordenadas Y', 'Fontsize', 12)
    % zticklabels('z = 0')
    grid on
    axis equal
    format shortE
    title(['Modo ' int2str(m) ' de{\it K},' ' autovalor ' num2str(autovalor, '%6.3e') ','...
    sprintf('\n') 'energía de deformación = ' num2str(PI_int(m), '%6.3e')], 'Fontsize', 12)
end
% saveas(gcf,'modos_deformacion.png')

%% Condiciones de contorno en desplazamientos (esenciales o de Dirichlet)
U = zeros(nGDL, 1); % vector inicial de desplazamientos
F = zeros(nGDL, 1);
gdlT = [1:nGDL]; % grados de libertad totales
gdlR = [1 2 3 nGDL-2 nGDL-1 nGDL]; % grados de libertad restringidos: SPC
uR = [0; 0; 0; 0; 0; 0]; % vector de desplazamientos impuestos: SPC
U(gdlR) = uR; % sustitución de desplazamientos impuestos con rot incluida
gdlL = setdiff(gdlT, gdlR);

%% Condensación estática
KLL = K(gdlL, gdlL);
KRR = K(gdlR, gdlR);
KLR = K(gdlL, gdlR);
KRL = KLR';
fL = FA(gdlL); % fuerzas externas en grados libres

% incógnitas
uL = KLL\(fL - KLR*uR)
fR = KRL*uL + KRR*uR % fuerzas gdlR en las restricciones

% Desplazamientos totales
U(gdlL) = uL;
Ux = U([1:3:21]);
Uy = U([2:3:21]);
Uz = U([3:3:21]);

% coordenadas de los nodos deformados
c_nod_def = c_nod + [Ux, Uy];

% Fuerzas totales
F(gdlL) = fL;
F(gdlR) = fR;
FR = fR - Fneq(gdlR) % fuerzas de reacción

%% Verificación
Fint = K*U;
Fext = F;
error_F = Fint - Fext; %verificación de equilibrio interno
if any(abs(error_F) > 1e-6) % Modificado ligeramente para prevenir error de matrices
    disp('error: no funciona')
else
    disp('Funciona OK')
end

%% Fuerzas y desplazamientos por elementos
U_e = zeros(6, nele);
F_e = zeros(6, nele);
Fint_e = zeros(6, nele); % Inicializar variable

for e = 1:nele
    gdl_e = conectividad_g(e,:);
    u_e = U(gdl_e); % desplazamientos en globales del ele1
    U_e(:,e) = u_e; % desplazamientos de los elementos
    Fint_e(:,e) = (K_ele(:,:,e)*u_e); % fuerzas nodales internas
end

%% Representación gráfica de la deformada
figure(1);
escala = L/max(abs([Ux; Uy]))*0.1; % escala de magnificación 10%
cnod_def_esc = c_nod + [Ux, Uy]*escala; % coordenadas de los nodos deformados escalados

% apaño rápido para elementos 1D
% plot(c_nod(:,1),c_nod(:,2),'b-o')
% hold on
% grid on
% axis equal
% title('Estructura original y deformada elástica','Fontsize',12)
% plot(cnoddef_esc(:,1),cnoddef_esc(:,2),'r--*')
% forma correcta para cualquier elemento

for e = 1:nele
    conecta_nodos = conectividad_n(e,:);
    % elemento sin deformar
    n1_d = cnod_def_esc(conecta_nodos(1),:); 
    n2_d = cnod_def_esc(conecta_nodos(2),:);
    plot([n1_d(1) n2_d(1)], [n1_d(2) n2_d(2)], 'r--o')
    % plot([n1(1)+u_e(1) n2(1)+u_e(4)],[n1(2)+u_e(2) n2(2)+u_e(5)],'b-*')
end
grid on
axis equal
title('Estructura original y deformada elástica', 'Fontsize', 12)

%% Representar el perfil o Skyline de la matriz de rigidez
% subplot(1,2,2)
figure(fn); fn = fn + 1;
verK = zeros(nGDL);

for s = 1:size(K,1)
    for t = 1:size(K,1)
        if abs(K(s,t)) > 1e-12
            verK(s,t) = 1;
        else 
            verK(s,t) = 0; % Corregido el original "verK(s,t)<1e-12;" que no asignaba valor
        end
    end
end
imagesc(verK); % Crea un plot de colores con los valores de K
colormap(flipud(gray)); % cambia el mapa de colores a gris
% los mayores negro, y los menores blanco
title('Ocupación de K', 'Fontsize', 10)
xlabel('gdln', 'Fontsize', 10) % Etiqueta el eje horizontal
ylabel('gdln', 'Fontsize', 10) % Etiqueta el eje vertical