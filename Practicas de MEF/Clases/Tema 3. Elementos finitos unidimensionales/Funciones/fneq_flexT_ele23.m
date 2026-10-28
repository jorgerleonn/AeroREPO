function [Fneq_f] = fneq_flexT_ele23(f_y,m_z,Le,c_i)
% Ke rigidez del elemento barra axil en coordenadas locales solo con
% grados de deformación, no están los de sólido rígido
% f_y vector columna con los valores de la función de carga en los nodos
% f_y = f0y*[1;1] % carga constante 2 nodos
% f_y = f0y*[0;1] % carga lineal 2 nodos
% f_y = f0y*[1;1;1] % carga constante 3 nodos
% f_y = f0y*[0;1/2;1] % carga lineal 3 nodos
% f_y = f0y*[0;1/4;1] % carga cuadrática 3 nodos
% nnodos número de nodos del elemento
% c_i cuadratura de integración de Ke
nnodos = length(f_y);
Je = Le/2; iJe = 1/Je;
if nnodos==2
    Fneq_f = zeros(4,1);
    [chi_ip,w_ip] = cuadratura(c_i);
    for i=1:length(chi_ip)
        chi = chi_ip(i); wi = w_ip(i);
        N1 = (1-chi)/2; N2 = (1+chi)/2;
        NL = [N1 N2]; % matriz cinemática del elemento 2nodos
        fy_chi = NL*f_y; % fuerza cortante interior del elemento virgulilla (aproximacion isoparam)
        mz_chi = NL*m_z; % momento flector interior del elemento virgulilla (aproximacion isoparam)
        Nf = [N1 0 N2 0]; % funciones de forma de cálculo de flecha
        Ng = [0 N1 0 N2]; % funciones de forma de cálculo de giro
        Fneq_f = Fneq_f + Nf'*fy_chi*Je*wi + Ng'*mz_chi*Je*wi; %este elemento es aquí defectuoso aquí: hay que corregirlo
        % "out of scope"
    end
elseif nnodos==3
    Fneq_f = zeros(6,1);
    [chi_ip,w_ip] = cuadratura(c_i);
    for i=1:length(chi_ip)
        chi = chi_ip(i); wi = w_ip(i);
        %
        Fneq_f = Fneq_f + Nf'*fy_chi*Je*wi + Ng'*mz_chi*Je*wi;
    end
else
    disp('error en el número de nodos')
    return
end
%
function [chi_ip,w_ip] = cuadratura(orden) % datos de integración numérica
if orden==1
    chi_ip = 0; w_ip = 2; %cuadratura n=1
elseif orden==2
    chi_ip = [-1 1]*1/sqrt(3); w_ip = [1 1]; %cuadratura n=2
elseif orden==3
    chi_ip = [-sqrt(3/5) 0 sqrt(3/5)]; w_ip = [5/9 8/9 5/9]; %cuadratura n=3
elseif orden==4
    chi_ip = [-sqrt((3+2*sqrt(6/5))/7) -sqrt((3-2*sqrt(6/5))/7) sqrt((3-2*sqrt(6/5))/7) sqrt((3+2*sqrt(6/5))/7)];
    w_ip = [(18-sqrt(30))/36 (18+sqrt(30))/36 (18+sqrt(30))/36 (18-sqrt(30))/36]; %cuadratura n=4
end
end
end