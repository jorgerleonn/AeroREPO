%% Problema 3 - Portico 3D con rotula en n5/n6, desplazamiento impuesto y carga distribuida
close all;

%% 1. Parametros
L = 1000;
E = 100e3;
nu = 1/3;
G = E/(2*(1+nu));

b = L/40;
a = 2*b;
t = L/400;
A = a*b - (a-2*t)*(b-2*t);
Iz = (b*a^3 - (b-2*t)*(a-2*t)^3)/12;
Iy = (a*b^3 - (a-2*t)*(b-2*t)^3)/12;
am = a - t;  bm = b - t;  Om = am*bm;
J = 4*Om^2 / (2*(am+bm)/t);

q0x    = 5;
u_load = L/20;

%% 2. Geometria (misma que P2; n6 duplicado de n5 -> rotula; resto +1)
c_nod = [ ...
   0,    0,   0;      % n1
   L/6,  0,   L/3;    % n2
   L/3,  0,   2*L/3;  % n3
   L/2,  0,   L;      % n4
   5*L/6,0,   L;      % n5
   5*L/6,0,   L;      % n6  (duplicado de n5 -> rotula)
   7*L/6,0,   L;      % n7
   3*L/2,0,   L;      % n8
   3*L/2,L/4, L;      % n9
   3*L/2,L/2, L;      % n10
   3*L/2,L/2, 2*L/3;  % n11
   3*L/2,L/2, L/3;    % n12
   3*L/2,L/2, 0;      % n13  (apoyo FR2)
   3*L/2,-L/4,L;      % n14
   3*L/2,-L/2,L;      % n15
   3*L/2,-L/2,2*L/3;  % n16
   3*L/2,-L/2,L/3;    % n17
   3*L/2,-L/2,0];     % n18  (apoyo)

conectividad_n = [ ...
   1  2;  2  3;  3  4;  4  5;      % columna + tramo hasta n5
   6  7;  7  8;  8  9;  9 10;      % n6 -> n7 -> n8 -> n9 -> n10
   8 14; 14 15;                    % rama -Y
   10 11; 11 12; 12 13;            % pata derecha
   15 16; 16 17; 17 18];           % pata izquierda

nnod = size(c_nod,1);
n_elementos = size(conectividad_n,1);
gdlxnod = 6;
n_grados_libertad = gdlxnod*nnod;

%% 3. Geometria de elementos y matrices de transformacion
vref_e = cell(n_elementos,1);
vref_e(1:3) = {[-1;0;0]};
vref_e(4:10) = {[0;0;1]};
vref_e(11:16) = {[0;-1;0]};

Le = zeros(n_elementos,1);
Te = zeros(12,12,n_elementos);
conectividad_g = zeros(n_elementos,12);
for i = 1:n_elementos
    nodos = conectividad_n(i,:);
    [Le(i), Te(:,:,i)] = Geom_viga3D(c_nod(nodos(1),:), c_nod(nodos(2),:), vref_e{i});
    conectividad_g(i,:) = [gdlxnod*(nodos(1)-1)+(1:gdlxnod), ...
                           gdlxnod*(nodos(2)-1)+(1:gdlxnod)];
end

%% 4. Matriz de rigidez global
K_quest = zeros(n_grados_libertad);
K_ele = zeros(12,12,n_elementos);
for i = 1:n_elementos
    gdl_i = conectividad_g(i,:);
    Ke = Ke_viga3D(E,G,Le(i),A,Iy,Iz,J);
    K_e = Te(:,:,i)'*Ke*Te(:,:,i);
    K_quest(gdl_i,gdl_i) = K_quest(gdl_i,gdl_i) + K_e;
    K_ele(:,:,i) = K_e;
end

%% 5. Carga distribuida: constante en el travesano (z=L) y triangular en las patas
Fneq_ele = zeros(12,n_elementos);
Fneq = zeros(n_grados_libertad,1);
for e = 1:n_elementos
    nA = conectividad_n(e,1);
    nB = conectividad_n(e,2);
    if ismember(e,[1 2 3 4 5 6])
        qA = 0; qB = 0;
    else
        qA = q0x*c_nod(nA,3)/L;
        qB = q0x*c_nod(nB,3)/L;
    end
    wA = [qA;0;0];
    wB = [qB;0;0];
    Re = Te(1:3,1:3,e);
    f_loc = cargas_dist_elem(Re,Le(e),wA,wB);
    Fneq_ele(:,e) = f_loc;
    gdl_e = conectividad_g(e,:);
    Fneq(gdl_e) = Fneq(gdl_e) + Te(:,:,e)'*f_loc;
end

%% 6. Condiciones de contorno, rotula (MPC) y desplazamiento impuesto
gdlT = 1:n_grados_libertad;

gdl_n1  = 6*(1-1)+(1:6);
gdl_n3  = 6*(3-1)+(1:6);
gdl_n5  = 6*(5-1)+(1:6);
gdl_n6  = 6*(6-1)+(1:6);
gdl_n13 = 6*(13-1)+(1:6);
gdl_n18 = 6*(18-1)+(1:6);

gdl_apoyo = [gdl_n1, gdl_n13, gdl_n18];

gdl_esclavo = gdl_n6(1:3);
gdl_maestro = gdl_n5(1:3);

gdl_impuesto = gdl_n3(1);

gdl_libres = setdiff(gdlT, [gdl_apoyo, gdl_esclavo, gdl_impuesto]);
n_libres = length(gdl_libres);

T = zeros(n_grados_libertad, n_libres);
T(gdl_libres,:) = eye(n_libres);
for k = 1:3
    T(gdl_esclavo(k),:) = T(gdl_maestro(k),:);
end

U0 = zeros(n_grados_libertad,1);
U0(gdl_impuesto) = u_load;

K_red = T'*K_quest*T;
F_red = T'*(Fneq - K_quest*U0);
q = K_red \ F_red;
U_quest = T*q + U0;

%% 7. Fuerzas de reaccion y fuerza en la rotula
R = K_quest*U_quest - Fneq;

FR1_quest = R(gdl_n1);
FR2_quest = R(gdl_n13);
F_rotula_quest = R(gdl_n6(1:3));

%% 8. Comprobacion
fprintf('residuo max en GDL libres = %.3e\n', max(abs(R(gdl_libres))));
fprintf('max|U| = %.4f mm\n', max(abs(U_quest)));
fprintf('u(n3) = %.4f  u(n5) = %s\n', U_quest(gdl_impuesto), mat2str(U_quest(gdl_n5(1:3)).',6));
fprintf('F_rotula(n6) = %s   F_rotula(n5) = %s\n', ...
        mat2str(F_rotula_quest.',6), mat2str(R(gdl_n5(1:3)).',6));

%% 9. Representacion de la estructura deformada
Ux = U_quest(1:6:end); Uy = U_quest(2:6:end); Uz = U_quest(3:6:end);
dim_est = max(max(c_nod)-min(c_nod));
escala = 0.07*dim_est/max(abs(U_quest));
c_def = [c_nod(:,1)+Ux*escala, c_nod(:,2)+Uy*escala, c_nod(:,3)+Uz*escala];

figure(1); hold on; grid on; axis equal;
for i = 1:n_elementos
    nodos = conectividad_n(i,:);
    plot3(c_nod(nodos,1), c_nod(nodos,2), c_nod(nodos,3), 'b-o', 'LineWidth', 1.5)
    plot3(c_def(nodos,1), c_def(nodos,2), c_def(nodos,3), 'r--o', 'LineWidth', 1.5)
end
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
title('Estructura sin deformar y deformada (escalada)');
view(3); hold off;

%% 10. Funciones
function [Ke] = Ke_viga3D(E,G,Le,A,Iy,Iz,J)
    k_ax = E*A/Le;  k_t = G*J/Le;
    az = 12*E*Iz/Le^3; bz = 6*E*Iz/Le^2; cz = 4*E*Iz/Le; dz = 2*E*Iz/Le;
    ay = 12*E*Iy/Le^3; by = 6*E*Iy/Le^2; cy = 4*E*Iy/Le; dy = 2*E*Iy/Le;
    Ke = zeros(12);
    Ke([1 7],[1 7]) = [k_ax -k_ax; -k_ax k_ax];
    Ke([4 10],[4 10]) = [k_t -k_t; -k_t k_t];
    Ke([2 6 8 12],[2 6 8 12]) = [az bz -az bz; bz cz -bz dz; -az -bz az -bz; bz dz -bz cz];
    Ke([3 5 9 11],[3 5 9 11]) = [ay -by -ay -by; -by cy by dy; -ay by ay by; -by dy by cy];
end

function f_loc = cargas_dist_elem(Re,Le,w1,w2)
    p = [Re*w1, Re*w2];
    p1 = p(:,1); dp = p(:,2)-p(:,1);
    fu1 = Le*(p1(1)/2 + dp(1)/6);  fu2 = Le*(p1(1)/2 + dp(1)/3);
    fv1 = Le*(p1(2)/2 + dp(2)*3/20); fv2 = Le*(p1(2)/2 + dp(2)*7/20);
    ftz1 = Le^2*(p1(2)/12 + dp(2)/30); ftz2 = Le^2*(-p1(2)/12 - dp(2)/20);
    fw1 = Le*(p1(3)/2 + dp(3)*3/20); fw2 = Le*(p1(3)/2 + dp(3)*7/20);
    fty1 = -Le^2*(p1(3)/12 + dp(3)/30); fty2 = Le^2*(p1(3)/12 + dp(3)/20);
    f_loc = [fu1; fv1; fw1; 0; fty1; ftz1; fu2; fv2; fw2; 0; fty2; ftz2];
end

function [Le,Te] = Geom_viga3D(n1,n2,v_ref)
    n1 = n1(:); n2 = n2(:);
    delta = n2-n1; Le = norm(delta); n_x = delta/Le;
    v_ref = v_ref(:);
    n_z = cross(n_x,v_ref); n_z = n_z/norm(n_z);
    n_y = cross(n_z,n_x);   n_y = n_y/norm(n_y);
    Re = [n_x'; n_y'; n_z'];
    R0 = zeros(3);
    Te = [Re R0 R0 R0; R0 Re R0 R0; R0 R0 Re R0; R0 R0 R0 Re];
end
