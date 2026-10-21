%% Grader 1
close all;

%% 1. Parámetros de la estructura y del material
L = 1000; % Longitud característica [mm]
E = 100e3; % módulo de young [N/mm^2]
b = L/40; % ancho de la sección de la viga [mm] (a = 2b = L/20)
a = 2*b; % altura de la sección de la viga [mm]
t = L/400; % espesor constante de la sección de la viga [mm]
[I, A] = inercia_viga_cajon_rectangular(a, b, t);
Iz = I(1);   % eje fuerte (a paralelo al eje local Y)
Iy = I(2);   % eje debil

nu = 1/3;                        % coeficiente de Poisson
G = E/(2*(1+nu));                % modulo de cortante
J = 2*t*(a-t)^2*(b-t)^2/(a+b-2*t);   % torsion, seccion cerrada de pared delgada (Bredt)

P1 = 10000; % Carga puntual en el nodo 8 [N]
P2 = 5000;  % Carga puntual en el nodo 13 [N]
q0x = 5; % Carga distribuida lineal [N/mm]

%% 2. generación de la geometría, malla de nodos y elementos, y grados de libertad
dirX = [1 0 0]'; dirY = [0 1 0]'; dirZ = [0 0 1]'; DIR = [dirX, dirY, dirZ];

% Coordenadas de los nodos [x, y, z]
n1  = [0, 0, 0];
n2  = [L/6, 0, L/3];
n3  = [L/3, 0, 2*L/3];
n4  = [L/2, 0, L];
n5  = [5*L/6, 0, L];
n6  = [7*L/6, 0, L];
n7  = [3*L/2, 0, L];
n8  = [3*L/2, L/4, L];
n9  = [3*L/2, L/2, L];
n10 = [3*L/2, L/2, 2*L/3];
n11 = [3*L/2, L/2, L/3];
n12 = [3*L/2, L/2, 0];
n13 = [3*L/2, -L/4, L];
n14 = [3*L/2, -L/2, L];
n15 = [3*L/2, -L/2, 2*L/3];
n16 = [3*L/2, -L/2, L/3];
n17 = [3*L/2, -L/2, 0];

c_nod = [n1; n2; n3; n4; n5; n6; n7; n8; n9; n10; n11; n12; n13; n14; n15; n16; n17]; % array de nodos
nnod = size(c_nod, 1); % número de nodos de la estructura

% Elementos y conectividad por nodos
conectividad_n = [1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 9; 13 7; 14 13; 10 9; 11 10; 12 11; 15 14; 16 15; 17 16]; % array de conectividad de los nodos (ejes locales x segun la figura: travesano en +Y y ramas verticales hacia arriba)
n_elementos = size(conectividad_n, 1); % número de elementos de la estructura

gdlxnod = 6; % Grados de libertad por nodo (3D)
n_nodos_elem = 2;
n_grados_libertad = gdlxnod * nnod; % número de grados de libertad

%% 2.1 Matrices de transformación y Longitud de los elementos
Le = zeros(n_elementos,1 ); % Longitud de los elementos
Te = zeros(12, 12, n_elementos);
conectividad_g = zeros(n_elementos, gdlxnod*n_nodos_elem); % Matriz de conectividad por grados de libertad

% Vector de referencia (eje local y del alma) de cada elemento segun el dibujo.
vref_e = cell(n_elementos, 1);
vref_e(1:3) = {[-1; 0; 0]}; % vigas inclinadas n1-n4 (y local = -X global)
vref_e(4:6) = {[0; 0; 1]}; % vigas rectas n4-n7 (y local = +Z global)
vref_e(7:8) = {[0; 0; 1]};  % rama +Y horizontal n7-n8-n9
vref_e(9:10) = {[0; 0; 1]}; % rama -Y horizontal n7-n13-n14
vref_e(11:13) = {[0; -1; 0]}; % rama +Y vertical n9-n12 (y local = -Y)
vref_e(14:16) = {[0; -1; 0]}; % rama -Y vertical n14-n17 (y local = -Y)


for i = 1:n_elementos
    nodos = conectividad_n(i,:);
    na = c_nod(nodos(1), :);
    nb = c_nod(nodos(2), :);

    % Longitud y matriz de transformacion del elemento (v_ref opcional)
    if isempty(vref_e{i})
        [Le(i), Te(:,:,i)] = Geom_viga3D(na, nb);
    else
        [Le(i), Te(:,:,i)] = Geom_viga3D(na, nb, vref_e{i});
    end

    % Conectividad por grados de libertad del elemento i
    gdl_i = [gdlxnod*(nodos(1)-1) + (1:gdlxnod), ...
             gdlxnod*(nodos(2)-1) + (1:gdlxnod)];
    conectividad_g(i, :) = gdl_i;
end

%% 3. matriz de rigidez global (consultar página 37 de las diapositias)
K_quest = zeros(n_grados_libertad, n_grados_libertad);
K_ele = zeros(12,12, n_elementos); % Matriz de rigidez de los elementos

for i = 1:n_elementos
% Ensamblaje de la matriz de rigidez global
    gdl_i = conectividad_g(i, :); % Grados de libertad del elemento i
    Te_i = Te(:,:,i); % Matriz de transformación del elemento i
    Ke_i = Ke_viga3D(E, G, Le(i), A, Iy, Iz, J); % Matriz de rigidez local del elemento i
    K_e = Te_i'*Ke_i*Te_i; % Matriz de rigidez global del elemento i
    K_quest(gdl_i, gdl_i) = K_quest(gdl_i, gdl_i) + K_e; % Ensamblaje de la matriz de rigidez global
    K_ele(:,:,i) = K_e; % Guardar matriz de rigidez del elemento i
end

%% 4. vector de fuerzas nodales equivalentes
Fneq_ele = zeros(12, n_elementos);
Fneq = zeros(n_grados_libertad, 1); % Vector de fuerzas nodales equivalentes globales

% Carga distribuida en +X global con perfil triangular q(z) = q0x * z / L:
%   - columna n1-n4: 0 en n1 y q0x en n4
%   - travesano superior n14-n9 (z = L): constante q0x
%   - ramas verticales n14-n17 y n9-n12: q0x arriba, 0 abajo
%   - la viga superior n4-n7 (elementos 4-6) NO lleva carga distribuida
for e = 1:n_elementos
    nA = conectividad_n(e, 1);
    nB = conectividad_n(e, 2);

    if ismember(e, [4 5 6])
        q_med = 0;
    else
        z_med = 0.5 * (c_nod(nA, 3) + c_nod(nB, 3)); % centroide del elemento segun Z
        q_med = q0x * z_med / L;                     % valor de q0x*z/L en el centroide
    end
    w_med = [q_med; 0; 0]; % carga uniforme (global) en el elemento e

    Re = Te(1:3, 1:3, e); % bloque de rotacion global -> local del elemento e
    f_loc = cargas_dist_elem(Re, Le(e), w_med, w_med);
    Fneq_ele(:, e) = f_loc;

    gdl_e = conectividad_g(e, :);
    Fneq(gdl_e) = Fneq(gdl_e) + Te(:,:,e)' * f_loc; % a ejes globales
end

% Cargas puntuales: P1 en -Z global en el nodo 8 y P2 en -Z global en el nodo 13
Fn = zeros(n_grados_libertad, 1); % Vector de fuerzas nodales puntuales globales
Fn(6*(8-1) + 3)  = -P1; % w (Z) del nodo 8
Fn(6*(13-1) + 3) = -P2; % w (Z) del nodo 13
F_total = Fn + Fneq; % Vector de fuerzas nodales equivalentes globales total (distribuidas + puntuales)

Fneq_quest = F_total; % fuerzas nodales equivalentes de la estructura (distribuidas + puntuales P1 y P2)

%% 5. imposición de las condiciones de contorno
U_quest = zeros(n_grados_libertad, 1); % Vector de desplazamientos nodales globales
gdlT = 1:n_grados_libertad; % Grados de libertad totales

% Apoyos empotrados (6 GDL) en n1 (columna), n12 y n17 (ramas)
gdl_apoyo = [1:6, 6*(12-1)+(1:6), 6*(17-1)+(1:6)];
gdlR = gdl_apoyo; % GDL con desplazamiento prescrito (apoyos)
gdl_libres = setdiff(gdlT, gdlR); % Grados de libertad libres

K_LL = K_quest(gdl_libres, gdl_libres); % Submatriz de rigidez global correspondiente a los grados de libertad libres
K_LR = K_quest(gdl_libres, gdlR); % Submatriz de rigidez global correspondiente a los grados de libertad libres y restringidos
K_RR = K_quest(gdlR, gdlR); % Submatriz de rigidez global correspondiente a los grados de libertad restringidos
K_RL = K_LR'; % Submatriz de rigidez global correspondiente a los grados de libertad restringidos y libres

u_R = zeros(length(gdlR), 1); % Desplazamientos restringidos (todos nulos)

%% 6. resolución de incógnitas del problema

% 6.1 desplazamientos de los nodos
u_l = K_LL \ (F_total(gdl_libres) - K_LR * u_R); % Resolver el sistema (incluye el efecto del desplazamiento impuesto)
U_quest(gdlR) = u_R; % Asignar desplazamientos restringidos al vector de desplazamientos nodales globales
U_quest(gdl_libres) = u_l; % Asignar desplazamientos libres al vector de desplazamientos nodales globales

Ux = U_quest(1:6:end); % Desplazamientos en x
Uy = U_quest(2:6:end); % Desplazamientos en y
Uz = U_quest(3:6:end); % Desplazamientos en z

% Coordenadas de los nodos deformados
c_nod_def = zeros(nnod, 2); % Inicializar matriz de coordenadas de los nodos deformados
for i = 1:nnod
    c_nod_def(i, 1) = c_nod(i, 1) + Ux(i); % Coordenada x del nodo deformado
    c_nod_def(i, 2) = c_nod(i, 2) + Uy(i); % Coordenada y del nodo deformado
    c_nod_def(i, 3) = c_nod(i, 3) + Uz(i); % Coordenada z del nodo deformado
end

% 6.2 fuerzas de reacción en SPCs
F_int = K_quest * U_quest; % Fuerzas internas globales de la estructura
FR1_quest  = F_int(1:6)    - F_total(1:6);    % Reacciones en el apoyo 1 (nodo 1)
FR2_quest  = F_int(97:102) - F_total(97:102); % Reacciones en el apoyo 2 (nodo 17)
FR12_quest = F_int(67:72)  - F_total(67:72);  % Reacciones en el apoyo n12 (auxiliar)

% 6.4 fuerzas en el elemento en coordenadas locales
F_ele = zeros(12, n_elementos); % Fuerzas nodales de cada elemento en ejes globales
Fe_loc = zeros(12, n_elementos); % Fuerzas nodales de cada elemento en ejes locales
for e = 1:n_elementos
    gdl_e = conectividad_g(e, :); % Grados de libertad del elemento e
    u_e = U_quest(gdl_e); % Desplazamientos del elemento en ejes globales
    F_ele(:, e) = K_ele(:,:,e) * u_e; % Fuerzas nodales del elemento (globales) debidas a la rigidez
    Fe_loc(:, e) = Te(:,:,e) * F_ele(:, e); % Fuerzas nodales del elemento en ejes locales
end

% Fuerzas internas locales de los elementos (columna [12x1])
Fe1_quest = Fe_loc(:, 1);
Fe4_quest = Fe_loc(:, 4);
Fe16_quest = Fe_loc(:, 16);

% Momento interno de torsión del elemento 10 (GDL local x, θx, nodo 1)
Mt10_quest = Fe_loc(4, 10);

%% 7. representaciones gráficas
% 7.1 Representación gráfica de la estructura sin deformar y deformada
figure(2);
hold on; grid on; axis equal;

% Factor de escala para que la deformada sea apreciable (15% de la dimensión)
dim_est = max(max(c_nod) - min(c_nod)); % dimensión característica de la estructura
escala = 0.10 * dim_est / max(abs(U_quest)); % evita que la deformada quede oculta
c_nod_def_esc = c_nod + [Ux, Uy, Uz] * escala; % nodos deformados escalados

for i = 1:n_elementos
    nodos = conectividad_n(i,:);
    na = c_nod(nodos(1),:); 
    nb = c_nod(nodos(2),:);
    na_def = c_nod_def_esc(nodos(1),:);
    nb_def = c_nod_def_esc(nodos(2),:);
    plot3([na(1) nb(1)], [na(2) nb(2)], [na(3) nb(3)], 'b-o', 'LineWidth', 1.5) % sin deformar
    plot3([na_def(1) nb_def(1)], [na_def(2) nb_def(2)], [na_def(3) nb_def(3)], 'r--o', 'LineWidth', 1.5) % deformada escalada
end

xlabel('X (m)','Fontsize', 11);
ylabel('Y (m)','Fontsize', 11);
zlabel('Z (m)','Fontsize', 11);
title('Estructura original y deformada elástica (deformada escalada)', 'Fontsize', 12)
legend('Sin deformar','Deformada','Location','best')
view(3)
hold off;

% 7.2 Representación de la ocupación de la matriz de rigidez global
figure(1);
verK = zeros(n_grados_libertad);
for i=1:size(K_quest,1)
    for j=1:size(K_quest,1)
        if abs(K_quest(i,j)) > 1e-12
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

%% 8. Subrutinas de usuario - functions

function [Ke] = Ke_viga3D(E, G, Le, A, Iy, Iz, J)
    % Rigidez local de una viga 3D de 2 nodos (12x12)
    % E  : modulo de Young
    % G  : modulo de cortante
    % Le : longitud del elemento
    % A  : area de la seccion
    % Iy : inercia de flexion respecto al eje local y
    % Iz : inercia de flexion respecto al eje local z
    % J  : constante de torsion
    % Orden de GDL: [u1 v1 w1 tx1 ty1 tz1 | u2 v2 w2 tx2 ty2 tz2]

    k_ax = E*A/Le;    % rigidez axial
    k_t  = G*J/Le;    % rigidez a torsion

    az = 12*E*Iz/Le^3; % flexion en plano x-y (v, tz)
    bz =  6*E*Iz/Le^2;
    cz =  4*E*Iz/Le;
    dz =  2*E*Iz/Le;

    ay = 12*E*Iy/Le^3; % flexion en plano x-z (w, ty)
    by =  6*E*Iy/Le^2;
    cy =  4*E*Iy/Le;
    dy =  2*E*Iy/Le;

    Ke = zeros(12, 12);

    % Axial (u1, u2)
    Ke([1 7], [1 7]) = [ k_ax -k_ax;
                        -k_ax  k_ax];

    % Torsion (tx1, tx2)
    Ke([4 10], [4 10]) = [ k_t -k_t;
                          -k_t  k_t];

    % Flexion en x-y: v1, tz1, v2, tz2
    Ke([2 6 8 12], [2 6 8 12]) = [ az  bz -az  bz;
                                   bz  cz -bz  dz;
                                  -az -bz  az -bz;
                                   bz  dz -bz  cz];

    % Flexion en x-z: w1, ty1, w2, ty2
    Ke([3 5 9 11], [3 5 9 11]) = [ ay -by -ay -by;
                                  -by  cy  by  dy;
                                  -ay  by  ay  by;
                                  -by  dy  by  cy];
end

function f_loc = cargas_dist_elem(Re, Le, w1, w2)
    % Fuerzas nodales equivalentes (consistentes) de una carga distribuida
    % lineal sobre un elemento viga 3D de 2 nodos.
    % w1, w2 : carga por unidad de longitud en los nodos 1 y 2, en EJES GLOBALES.
    % Re     : matriz 3x3 de rotacion global -> local del elemento.
    p = [Re*w1, Re*w2];      % carga en ejes locales en nodo 1 y nodo 2
    p1 = p(:,1);
    dp = p(:,2) - p(:,1);

    % Axial (u1, u2)
    fu1 = Le*(p1(1)/2 + dp(1)/6);
    fu2 = Le*(p1(1)/2 + dp(1)/3);

    % Flexion en el plano local x-y (v, tz)
    fv1  = Le*(p1(2)/2 + dp(2)*3/20);
    fv2  = Le*(p1(2)/2 + dp(2)*7/20);
    ftz1 = Le^2*( p1(2)/12 + dp(2)/30);
    ftz2 = Le^2*(-p1(2)/12 - dp(2)/20);

    % Flexion en el plano local x-z (w, ty)
    fw1  = Le*(p1(3)/2 + dp(3)*3/20);
    fw2  = Le*(p1(3)/2 + dp(3)*7/20);
    fty1 = -Le^2*( p1(3)/12 + dp(3)/30);
    fty2 =  Le^2*( p1(3)/12 + dp(3)/20);

    f_loc = [fu1; fv1; fw1; 0; fty1; ftz1; ...
             fu2; fv2; fw2; 0; fty2; ftz2];
end

function [Le, Te] = Geom_viga3D(n1, n2, v_ref)
    % n1, n2: coordenadas globales de los nodos locales 1 y 2
    % v_ref : vector (en globales) que orienta la sección transversal.
    %         Marca la direccion aproximada del eje local y (el "arriba"
    %         de la seccion). Si se omite, se toma [0;0;1] salvo en barras
    %         verticales, donde Z global no sirve.
    n1 = n1(:); n2 = n2(:);

    % 1) Longitud del elemento
    delta = n2 - n1;
    Le = norm(delta);

    % 2) Eje local x: a lo largo de la barra
    n_x = delta / Le;

    % 3) Vector de referencia por defecto
    if nargin < 3 || isempty(v_ref)
        if abs(n_x(3)) > 0.999    % barra vertical -> Z global paralelo a x
            v_ref = [0; -1; 0];
        else
            v_ref = [0; 0; 1];
        end
    else
        v_ref = v_ref(:);
    end

    % 4) Ejes locales z (perpendicular a x y a v_ref) e y (dextrogiro)
    n_z = cross(n_x, v_ref);
    if norm(n_z) < 1e-12
        error('Geom_viga3D: v_ref es paralelo al eje de la barra.');
    end
    n_z = n_z / norm(n_z);
    n_y = cross(n_z, n_x);
    n_y = n_y / norm(n_y);        % x e z ortonormales -> y ya unitario

    % 5) Matriz de cambio de base global -> local.
    %    Base global = identidad, luego las filas de Re son los ejes
    %    locales expresados en globales.
    Re = [n_x'; n_y'; n_z'];

    % 6) Transformacion 12x12 = blkdiag(Re,Re,Re,Re)
    R0 = zeros(3);
    Te = [Re R0 R0 R0;
          R0 Re R0 R0;
          R0 R0 Re R0;
          R0 R0 R0 Re];
end

function [I, A] = inercia_viga_cajon_rectangular(a, b, t)
    % a: altura total de la viga
    % b: Anchura total de las alas (base)
    % t: espesor constante de las almas y las alas
    % I: [Iy, Iz] momentos de inercia de la seccion (Iy = eje fuerte)
    % A: area de la seccion
    Iy = (b * a^3 - (b - 2*t) * (a - 2*t)^3) / 12;   % eje fuerte (a segun y)
    Iz = (a * b^3 - (a - 2*t) * (b - 2*t)^3) / 12;   % eje debil
    A = a * b - (a - 2*t) * (b - 2*t);

    I = [Iy, Iz];
end

% Fin del código