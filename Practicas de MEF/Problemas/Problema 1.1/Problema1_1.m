%% Grader 1

L = 1000;
E = 100e3;
b = L/40;
a = 2*b;
t = L/400;
[I, A] = inercia_viga_cajon_rectangular(a, b, t);
Iz = I(1);
Iy = I(2);

nu = 1/3;
G = E/(2*(1+nu));
J = 2*t*(a-t)^2*(b-t)^2/(a+b-2*t);

P1 = 10000;
P2 = 5000;
q0x = 5;

%% 2. generación de la geometría, malla de nodos y elementos, y grados de libertad
dirX = [1 0 0]'; dirY = [0 1 0]'; dirZ = [0 0 1]'; DIR = [dirX, dirY, dirZ];

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

c_nod = [n1; n2; n3; n4; n5; n6; n7; n8; n9; n10; n11; n12; n13; n14; n15; n16; n17];
nnod = size(c_nod, 1);

conectividad_n = [1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 9; 13 7; 14 13; 10 9; 11 10; 12 11; 15 14; 16 15; 17 16];
n_elementos = size(conectividad_n, 1);

gdlxnod = 6;
n_nodos_elem = 2;
n_grados_libertad = gdlxnod * nnod;

%% 2.1 Matrices de transformación y Longitud de los elementos
Le = zeros(n_elementos,1 );
Te = zeros(12, 12, n_elementos);
conectividad_g = zeros(n_elementos, gdlxnod*n_nodos_elem);

vref_e = cell(n_elementos, 1);
vref_e(1:3) = {[-1; 0; 0]};
vref_e(4:6) = {[0; 0; 1]};
vref_e(7:8) = {[0; 0; 1]};
vref_e(9:10) = {[0; 0; 1]};
vref_e(11:13) = {[0; -1; 0]};
vref_e(14:16) = {[0; -1; 0]};

for i = 1:n_elementos
    nodos = conectividad_n(i,:);
    na = c_nod(nodos(1), :);
    nb = c_nod(nodos(2), :);

    if isempty(vref_e{i})
        [Le(i), Te(:,:,i)] = Geom_viga3D(na, nb);
    else
        [Le(i), Te(:,:,i)] = Geom_viga3D(na, nb, vref_e{i});
    end

    gdl_i = [gdlxnod*(nodos(1)-1) + (1:gdlxnod), ...
             gdlxnod*(nodos(2)-1) + (1:gdlxnod)];
    conectividad_g(i, :) = gdl_i;
end

%% 3. matriz de rigidez global
K_quest = zeros(n_grados_libertad, n_grados_libertad);
K_ele = zeros(12,12, n_elementos);

for i = 1:n_elementos

    gdl_i = conectividad_g(i, :);
    Te_i = Te(:,:,i);
    Ke_i = Ke_viga3D(E, G, Le(i), A, Iy, Iz, J);
    K_e = Te_i'*Ke_i*Te_i;
    K_quest(gdl_i, gdl_i) = K_quest(gdl_i, gdl_i) + K_e;
    K_ele(:,:,i) = K_e;
end

%% 4. vector de fuerzas nodales equivalentes
Fneq_ele = zeros(12, n_elementos);
Fneq = zeros(n_grados_libertad, 1);

for e = 1:n_elementos
    nA = conectividad_n(e, 1);
    nB = conectividad_n(e, 2);

    if ismember(e, [4 5 6])
        qA = 0; qB = 0;
    else
        qA = q0x * c_nod(nA,3)/L;
        qB = q0x * c_nod(nB,3)/L;
    end
    wA = [qA; 0; 0];
    wB = [qB; 0; 0];

    Re = Te(1:3, 1:3, e);
    f_loc = cargas_dist_elem(Re, Le(e), wA, wB);
    Fneq_ele(:, e) = f_loc;

    gdl_e = conectividad_g(e, :);
    Fneq(gdl_e) = Fneq(gdl_e) + Te(:,:,e)' * f_loc;
end

Fn = zeros(n_grados_libertad, 1);
Fn(6*(8-1) + 3)  = -P1;
Fn(6*(13-1) + 3) = -P2;
F_total = Fn + Fneq;

Fneq_quest = Fneq;

%% 5. imposición de las condiciones de contorno
U_quest = zeros(n_grados_libertad, 1);
gdlT = 1:n_grados_libertad;

gdl_apoyo = [1:6, 6*(12-1)+(1:6), 6*(17-1)+(1:6)];
gdlR = gdl_apoyo;
gdl_libres = setdiff(gdlT, gdlR);

K_LL = K_quest(gdl_libres, gdl_libres);
K_LR = K_quest(gdl_libres, gdlR);
K_RR = K_quest(gdlR, gdlR);
K_RL = K_LR';

u_R = zeros(length(gdlR), 1);

%% 6. resolución de incógnitas del problema

u_l = K_LL \ (F_total(gdl_libres) - K_LR * u_R);
U_quest(gdlR) = u_R;
U_quest(gdl_libres) = u_l;

Ux = U_quest(1:6:end);
Uy = U_quest(2:6:end);
Uz = U_quest(3:6:end);

c_nod_def = zeros(nnod, 2);
for i = 1:nnod
    c_nod_def(i, 1) = c_nod(i, 1) + Ux(i);
    c_nod_def(i, 2) = c_nod(i, 2) + Uy(i);
    c_nod_def(i, 3) = c_nod(i, 3) + Uz(i);
end

F_int = K_quest * U_quest;
FR1_quest = F_int(1:6)    - F_total(1:6);
FR2_quest = F_int(67:72)  - F_total(67:72);

F_ele = zeros(12, n_elementos);
Fe_loc = zeros(12, n_elementos);
for e = 1:n_elementos
    gdl_e = conectividad_g(e, :);
    u_e = U_quest(gdl_e);
    F_ele(:, e) = K_ele(:,:,e) * u_e;
    Fe_loc(:, e) = Te(:,:,e) * F_ele(:, e);
end

Fe1_quest = Fe_loc(:, 1);
Fe4_quest = Fe_loc(:, 4);
Fe16_quest = Fe_loc(:, 16);

Mt10_quest = Fe_loc(4, 10);

%% 7. representaciones gráficas

figure(2);
hold on; grid on; axis equal;

dim_est = max(max(c_nod) - min(c_nod));
escala = 0.10 * dim_est / max(abs(U_quest));
c_nod_def_esc = c_nod + [Ux, Uy, Uz] * escala;

for i = 1:n_elementos
    nodos = conectividad_n(i,:);
    na = c_nod(nodos(1),:);
    nb = c_nod(nodos(2),:);
    na_def = c_nod_def_esc(nodos(1),:);
    nb_def = c_nod_def_esc(nodos(2),:);
    plot3([na(1) nb(1)], [na(2) nb(2)], [na(3) nb(3)], 'b-o', 'LineWidth', 1.5)
    plot3([na_def(1) nb_def(1)], [na_def(2) nb_def(2)], [na_def(3) nb_def(3)], 'r--o', 'LineWidth', 1.5)
end

xlabel('X (m)','Fontsize', 11);
ylabel('Y (m)','Fontsize', 11);
zlabel('Z (m)','Fontsize', 11);
title('Estructura original y deformada elástica (deformada escalada)', 'Fontsize', 12)
legend('Sin deformar','Deformada','Location','best')
view(3)
hold off;

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

    k_ax = E*A/Le;
    k_t  = G*J/Le;

    az = 12*E*Iz/Le^3;
    bz =  6*E*Iz/Le^2;
    cz =  4*E*Iz/Le;
    dz =  2*E*Iz/Le;

    ay = 12*E*Iy/Le^3;
    by =  6*E*Iy/Le^2;
    cy =  4*E*Iy/Le;
    dy =  2*E*Iy/Le;

    Ke = zeros(12, 12);

    Ke([1 7], [1 7]) = [ k_ax -k_ax;
                        -k_ax  k_ax];

    Ke([4 10], [4 10]) = [ k_t -k_t;
                          -k_t  k_t];

    Ke([2 6 8 12], [2 6 8 12]) = [ az  bz -az  bz;
                                   bz  cz -bz  dz;
                                  -az -bz  az -bz;
                                   bz  dz -bz  cz];

    Ke([3 5 9 11], [3 5 9 11]) = [ ay -by -ay -by;
                                  -by  cy  by  dy;
                                  -ay  by  ay  by;
                                  -by  dy  by  cy];
end

function f_loc = cargas_dist_elem(Re, Le, w1, w2)

    p = [Re*w1, Re*w2];
    p1 = p(:,1);
    dp = p(:,2) - p(:,1);

    fu1 = Le*(p1(1)/2 + dp(1)/6);
    fu2 = Le*(p1(1)/2 + dp(1)/3);

    fv1  = Le*(p1(2)/2 + dp(2)*3/20);
    fv2  = Le*(p1(2)/2 + dp(2)*7/20);
    ftz1 = Le^2*( p1(2)/12 + dp(2)/30);
    ftz2 = Le^2*(-p1(2)/12 - dp(2)/20);

    fw1  = Le*(p1(3)/2 + dp(3)*3/20);
    fw2  = Le*(p1(3)/2 + dp(3)*7/20);
    fty1 = -Le^2*( p1(3)/12 + dp(3)/30);
    fty2 =  Le^2*( p1(3)/12 + dp(3)/20);

    f_loc = [fu1; fv1; fw1; 0; fty1; ftz1; ...
             fu2; fv2; fw2; 0; fty2; ftz2];
end

function [Le, Te] = Geom_viga3D(n1, n2, v_ref)

    n1 = n1(:); n2 = n2(:);

    delta = n2 - n1;
    Le = norm(delta);

    n_x = delta / Le;

    if nargin < 3 || isempty(v_ref)
        if abs(n_x(3)) > 0.999
            v_ref = [0; -1; 0];
        else
            v_ref = [0; 0; 1];
        end
    else
        v_ref = v_ref(:);
    end

    n_z = cross(n_x, v_ref);
    if norm(n_z) < 1e-12
        error('Geom_viga3D: v_ref es paralelo al eje de la barra.');
    end
    n_z = n_z / norm(n_z);
    n_y = cross(n_z, n_x);
    n_y = n_y / norm(n_y);

    Re = [n_x'; n_y'; n_z'];

    R0 = zeros(3);
    Te = [Re R0 R0 R0;
          R0 Re R0 R0;
          R0 R0 Re R0;
          R0 R0 R0 Re];
end

function [I, A] = inercia_viga_cajon_rectangular(a, b, t)

    Iy = (b * a^3 - (b - 2*t) * (a - 2*t)^3) / 12;
    Iz = (a * b^3 - (a - 2*t) * (b - 2*t)^3) / 12;
    A = a * b - (a - 2*t) * (b - 2*t);

    I = [Iy, Iz];
end
