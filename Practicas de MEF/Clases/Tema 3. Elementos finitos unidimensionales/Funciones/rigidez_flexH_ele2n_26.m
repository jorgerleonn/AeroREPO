function [Kef] = rigidez_flexH_ele2n_26(Le, I_x, E, GL)
% RIGIDEZ_FLEXH_ELE2N_26 Matriz de rigidez a flexión de viga Euler-Bernoulli (2 nodos)
% Emplea polinomios cúbicos de Hermite y cuadratura de Gauss-Legendre 1D.
%
% Grados de libertad locales: [w1, theta1, w2, theta2]'
%
% Sintaxis:
%   [Kef] = rigidez_flexH_ele2n_26(Le, I_x, E, GL)
%
% Entradas:
%   Le   - Longitud del elemento [m o mm]
%   I_x  - Vector con el momento de inercia en los 2 nodos: [I1; I2]
%   E    - Módulo de elasticidad longitudinal [Pa o MPa]
%   GL   - Orden de cuadratura de Gauss-Legendre (1 a 5)
%
% Salida:
%   Kef  - Matriz de rigidez elemental (4x4)

    I_x = I_x(:); % Asegurar vector columna
    if length(I_x) ~= 2
        error('El vector de inercia I_x debe contener exactamente 2 nodos: [I1; I2].');
    end

    % Jacobiano 1D constante: dx/dchi = Le / 2
    Je = Le / 2;
    iJe2 = (1 / Je)^2;

    Kef = zeros(4);
    [chi_ip, w_ip] = cuadratura(GL);

    for ip = 1:length(chi_ip)
        chi = chi_ip(ip);
        wi  = w_ip(ip);

        % Derivadas segundas de las funciones cúbicas de Hermite respecto a chi: d2N/dchi2
        Bf1 = (3/2) * chi;
        Bf2 = (3/2) * chi - 1/2;
        Bf3 = -(3/2) * chi;
        Bf4 = (3/2) * chi + 1/2;

        % Matriz cinemática de curvatura: kappa = Bef * d_e
        % donde theta_i = dw/dx, por lo que las columnas de giro van escaladas por Je
        Bef = iJe2 * [Bf1, Je * Bf2, Bf3, Je * Bf4];

        % Interpolación lineal isoparamétrica de la inercia en el elemento
        Ne = [(1 - chi) / 2, (1 + chi) / 2];
        I_chi = Ne * I_x;

        % Integración numérica: int(Bef' * E * I * Bef * det(J) * dchi)
        Kef = Kef + (Bef' * E * I_chi * Bef) * (Je * wi);
    end
end

%% Subrutina local: Cuadratura de Gauss-Legendre 1D en [-1, 1]
function [chi_ip, w_ip] = cuadratura(orden)
    switch orden
        case 1
            chi_ip = 0; 
            w_ip   = 2;
            
        case 2
            p      = 1 / sqrt(3);
            chi_ip = [-p, p]; 
            w_ip   = [1, 1];
            
        case 3
            p      = sqrt(3/5);
            chi_ip = [-p, 0, p]; 
            w_ip   = [5/9, 8/9, 5/9];
            
        case 4
            p1     = sqrt((3 - 2*sqrt(6/5)) / 7);
            p2     = sqrt((3 + 2*sqrt(6/5)) / 7);
            w1     = (18 + sqrt(30)) / 36;
            w2     = (18 - sqrt(30)) / 36;
            chi_ip = [-p2, -p1, p1, p2];
            w_ip   = [w2, w1, w1, w2];
            
        case 5
            p1     = (1/3) * sqrt(5 - 2*sqrt(10/7));
            p2     = (1/3) * sqrt(5 + 2*sqrt(10/7));
            w1     = (322 + 13*sqrt(70)) / 900;
            w2     = (322 - 13*sqrt(70)) / 900;
            chi_ip = [-p2, -p1, 0, p1, p2];
            w_ip   = [w2, w1, 128/225, w1, w2];
            
        otherwise
            error('Orden de cuadratura no soportado. Usa un valor entre 1 y 5.');
    end
end