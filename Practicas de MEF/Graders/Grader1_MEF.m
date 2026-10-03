%% Grader 1
close all;

%% 1. Parámetros de la estructura y del material
L = 1; % Longitud característica
E = 100e9; % módulo de young en GPa
b = L/40; % ancho de la sección de la viga
a = 2*b; % altura de la sección de la viga
t = L/200; % espesor constante de la sección de la viga
[I, A] = inercia_viga_doble_T(a, b, t);
Ix = I(1);   % eje fuerte (flexion en el plano x-y local)
Iy = I(2);   % eje debil  (flexion en el plano x-z local)

nu = 0.3;                        % coeficiente de Poisson
G = E/(2*(1+nu));                % modulo de cortante
J = (2*b*t^3 + (a-2*t)*t^3)/3;   % torsion, seccion abierta de pared delgada
Iz = Ix;                         % eje fuerte -> eje local z

P = 1000; % Carga puntual en el nodo 5, 1000 Newtons
q0y = 500; % Carga distribuida lineal de valor 500 N/m
uload = L/100; % Desplazamiento impuesto
gdl_uload = 20; % u=19, v=20, w=21, θx=22, θy=23, θz=24 -> grados de libertad del nodo 4

%% 2. generación de la geometría, malla de nodos y elementos, y grados de libertad
dirX = [1 0 0]'; dirY = [0 1 0]'; dirZ = [0 0 1]'; DIR = [dirX, dirY, dirZ];

L1 = 1*L; L2 = 2*L; L3 = 4*L; L4 = 3*L;
n1 = [L1, 0, 0];
n2 = [L1, 0, L3/3];
n3 = [L1, 0, 2*(L3/3)];
n4 = [L1, 0 ,L3];
n5 = [L1, L2/2, (L3-L4)/2 + L4];
n6 = [L1, L2, L4];
n7 = [0, L2, L4];

c_nod = [n1; n2; n3; n4; n5; n6; n7]; % array de nodos
nnod = size(c_nod, 1); % número de nodos de la estructura

% Elementos y conectividad por nodos
conectividad_n = [1 2; 2 3; 3 4; 4 5; 5 6; 6 7]; % array de conectividad de los nodos
n_elementos = size(conectividad_n, 1); % número de elementos de la estructura

gdlxnod = 6; % Grados de libertad por nodo (3D)
n_nodos_elem = 2;
n_grados_libertad = gdlxnod * nnod; % número de grados de libertad

%% 2.1 Matrices de transformación y Longitud de los elementos
Le = zeros(n_elementos,1 ); % Longitud de los elementos
Te = zeros(12, 12, n_elementos);
conectividad_g = zeros(n_elementos, gdlxnod*n_nodos_elem); % Matriz de conectividad por grados de libertad

% Vector de referencia (eje local y del alma) de cada elemento segun el dibujo.
% Dejar vacio [] para usar el valor por defecto de Geom_viga3D.
vref_e = cell(n_elementos, 1);
vref_e(1:3) = {[0; -1; 0]};   % columna vertical n1-n4

comp_alpha = atan((L3-L4)/L2);
vref_e{4} = [0; sin(comp_alpha); cos(comp_alpha)];
vref_e{5} = vref_e{4};

vref_e{6} = [0; 0; 1];

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

%% 4. vector de fuerzas nodales equivalentes (consultar página 30 de las diapositias)
Fneq_ele = zeros(12, n_elementos);
Fneq = zeros(n_grados_libertad, 1); % Vector de fuerzas nodales equivalentes globales

% Cálculo de las Fuerzas Nodales Equivalentes (cortantes y momentos) para los tres elementos 
% que conforman la columna vertical del pórtico sometidos a la carga distribuida lateral.
for i = 1:3
    nodos = conectividad_n(i, :);
    q1 = q0y * c_nod(nodos(1),3)/L3; % intensidad global (+Y), positiva
    q2 = q0y * c_nod(nodos(2),3)/L3;

    Le_i = Le(i);
    dq = q2 - q1; % parte triangular

    % Carga distribuida = uniforme q1 + triangular dq (formulas estandar)
    V1 = q1*Le_i/2 + dq*3*Le_i/20;
    M1 = q1*Le_i^2/12 + dq*Le_i^2/30;
    V2 = q1*Le_i/2 + dq*7*Le_i/20;
    M2 = (q1*Le_i^2/12 + dq*Le_i^2/20);

    % +Y global = -y local -> signo global del vector
    f_loc = [0; -V1; 0; 0; 0; -M1; 0; -V2; 0; 0; 0; M2];
    Fneq_ele(:, i) = f_loc;

    gdl_i = conectividad_g(i, :);
    Fneq(gdl_i) = Fneq(gdl_i) + Te(:,:,i)' * f_loc; % Vector de cargas nodales equivalentes en ejes globales
end

Fneq_quest = Fneq';

Fn = zeros(n_grados_libertad, 1); % Vector de fuerzas nodales equivalentes globales
Fn(27) = P; % carga aplicada en el eje +Z global en el nodo 5
F_total = Fn + Fneq; % Vector de fuerzas nodales equivalentes globales total; 

%% 5. imposición de las condiciones de contorno
U_quest = zeros(n_grados_libertad, 1); % Vector de desplazamientos nodales globales
gdlT = 1:n_grados_libertad; % Grados de libertad totales
gdl_apoyo = [1, 2, 3, ... 
        n_grados_libertad-5, n_grados_libertad-4, n_grados_libertad-3, ...
        n_grados_libertad-2, n_grados_libertad-1, n_grados_libertad]; % Grados de libertad de los apoyos
gdlR = [gdl_apoyo, gdl_uload]; % GDL con desplazamiento prescrito (apoyos + uload)
gdl_libres = setdiff(gdlT, gdlR); % Grados de libertad libres

K_LL = K_quest(gdl_libres, gdl_libres); % Submatriz de rigidez global correspondiente a los grados de libertad libres
K_LR = K_quest(gdl_libres, gdlR); % Submatriz de rigidez global correspondiente a los grados de libertad libres y restringidos
K_RR = K_quest(gdlR, gdlR); % Submatriz de rigidez global correspondiente a los grados de libertad restringidos
K_RL = K_LR'; % Submatriz de rigidez global correspondiente a los grados de libertad restringidos y libres

u_R = zeros(length(gdlR), 1); % Desplazamientos restringidos
u_R(gdlR == gdl_uload) = uload;

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
FR1_quest = F_int(1:3) - F_total(1:3); % Reacciones en el apoyo 1 (GDL 1,2,3), columna
FR7_quest = F_int(37:42) - F_total(37:42); % Reacciones en el apoyo 7 (GDL 37..42), columna
% 6.3 fuerza equivalente del desplazamiento impuesto uload
F_uload = F_int(gdl_uload) - F_total(gdl_uload); % Reacción en el GDL impuesto, escalar [1x1]

% 6.4 fuerzas en el elemento en coordenadas locales
F_ele = zeros(12, n_elementos); % Fuerzas nodales de cada elemento en ejes globales
Fe_loc = zeros(12, n_elementos); % Fuerzas nodales de cada elemento en ejes locales
for e = 1:n_elementos
    gdl_e = conectividad_g(e, :); % Grados de libertad del elemento e
    u_e = U_quest(gdl_e); % Desplazamientos del elemento en ejes globales
    F_ele(:, e) = K_ele(:,:,e) * u_e; % Fuerzas nodales totales del elemento (globales)
    Fe_loc(:, e) = Te(:,:,e) * F_ele(:, e); % Fuerzas nodales del elemento en ejes locales
end

% Fuerzas internas locales de los elementos como vectores fila
Fe1_quest = Fe_loc(:, 1)';
Fe4_quest = Fe_loc(:, 4)';
Fe6_quest = Fe_loc(:, 6)';

%% 7. representaciones gráficas
% 7.1 Representación gráfica de la estructura sin deformar y deformada
figure(1);
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
figure(3);
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

function [x, iterations] = GaussSeidel(A, b, tol, itmax)
    n = length(b);
    x = zeros(n, 1);
    normVal = 1e6;
    
    % * _*Tolerence for method*_
    if nargin == 2
        tol = 1e-8; 
        itmax = 2000;
    end
    
    iter = 0;
    while normVal > tol && iter < itmax
        x_old = x;
        for i = 1:n
            sigma = 0;
            for j = 1:i-1
                sigma = sigma + A(i,j)*x(j);
            end
            for j = i+1:n
                sigma = sigma + A(i,j)*x_old(j);
            end
            x(i) = (1/A(i,i))*(b(i) - sigma);
        end
        iter = iter + 1;
        normVal = norm(x_old - x);
    end
    iterations = iter;
end

function [I, A] = inercia_viga_doble_T(a, b, t)
    % a: altura total de la viga
    % b: Anchura total de las alas (base)
    % t: espesor constante de las almas y las alas
    % I: [Ix, Iy] momentos de inercia de la seccion
    % A: area de la seccion

    Ix = (b*a^3 - (b-t)*(a-2*t)^3) / 12;
    Iy = (2*(t*b^3) + (a - 2*t)*t^3) / 12;

    % area = 2 alas + alma
    A = 2*b*t + (a - 2*t)*t;

    I = [Ix, Iy];
end

% Fin del código