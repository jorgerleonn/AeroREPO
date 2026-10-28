function [Fneq_a] = fneq_axil_ele23(f_x,Le,c_i)
% Ke rigidez del elemento barra axil en coordenadas locales solo con
% grados de deformación, no están los de sólido rígido
% f_x vector columna con los valores de la función de carga en los nodos
% f_x = f0x*[1;1] % carga constante 2 nodos
% f_x = f0x*[0;1] % carga lineal 2 nodos
% f_x = f0x*[1;1;1] % carga constante 3 nodos
% f_x = f0x*[0;1/2;1] % carga lineal 3 nodos
% f_x = f0x*[0;1/4;1] % carga cuadrática 3 nodos
% nnodos número de nodos del elemento
% c_i cuadratura de integración de Ke
nnodos = length(f_x);
Je = Le/2; iJe = 1/Je;
if nnodos==2
    Fneq_a = zeros(2,1);
    [chi_ip,w_ip] = cuadratura(c_i);
    for i=1:length(chi_ip)
        chi = chi_ip(i);
        wi = w_ip(i); Ne1 = (1-chi)/2; Ne2 = (1+chi)/2;
        Ne = [Ne1 Ne2]; % funciones de forma del elemento de 2 nodos
        f_chi = Ne*f_x; % aproximación isoparamética de fuerza aplicada
        Fneq_a = Fneq_a + Ne'*f_chi*Je*wi; %sumatorio
    end
elseif nnodos==3
    Fneq_a = zeros(3,1);
    [chi_ip,w_ip] = cuadratura(c_i);
    for i=1:length(chi_ip)
        chi = chi_ip(i);
        wi = w_ip(i);
        Ne1 = chi*(chi-1)/2; Ne2 = 1-chi^2; Ne3 = chi*(chi+1)/2;
        Ne = [Ne1 Ne2 Ne3]; % funciones de forma del elemento de 3 nodos
        f_chi = Ne*f_x; % aproximación isoparamética de fuerza aplicada
        Fneq_a = Fneq_a + Ne'*f_chi*Je*wi; %sumatorio
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
```[cite: 3]