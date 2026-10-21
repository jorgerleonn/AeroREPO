function [Kea] = rigidez_axil_e33n_26(Le, A_x, E, GL)
% RIGIDEZ_AXIL_E33N_26 Matriz de rigidez axil elemental (2 o 3 nodos)
%
% Sintaxis:
%   [Kea] = rigidez_axil_e33n_26(Le, A_x, E, GL)
%
% Entradas:
%   Le   - Longitud del elemento [mm o m]
%   A_x  - Vector columna o fila con el área en los nodos:
%          Elemento lineal (2 nodos):   A_x = [A1; A2]
%          Elemento cuadrático (3 nodos): A_x = [A1; A2; A3]
%   E    - Módulo de Young del material [MPa o Pa]
%   GL   - Orden de la cuadratura de Gauss-Legendre (1 a 5)
%
% Salida:
%   Kea  - Matriz de rigidez elemental (2x2 o 3x3) en coordenadas locales

    A_x = A_x(:); % Asegurar vector columna
    num_nodos = length(A_x);
    
    % Jacobiano 1D constante: dx/dchi = Le / 2
    Je = Le / 2;
    iJe = 1 / Je;
    
    [chi_ip, w_ip] = cuadratura(GL);
    Kea = zeros(num_nodos);

    switch num_nodos
        case 2
            % Elemento lineal de 2 nodos: chi in [-1, 1]
            % dN/dchi constante
            Be = iJe * [-0.5, 0.5];
            
            for i = 1:length(chi_ip)
                chi = chi_ip(i);
                wi  = w_ip(i);
                
                Ne = [(1 - chi)/2, (1 + chi)/2];
                A_chi = Ne * A_x; % Interpolación isoparamétrica del área
                
                Kea = Kea + (Be' * E * A_chi * Be) * (Je * wi);
            end

        case 3
            % Elemento cuadrático de 3 nodos: chi in [-1, 1] (nodo central chi = 0)
            for i = 1:length(chi_ip)
                chi = chi_ip(i);
                wi  = w_ip(i);
                
                % Derivadas dN/dx = (1/Je) * dN/dchi
                Be = iJe * [chi - 0.5, -2 * chi, chi + 0.5];
                
                Ne = [chi * (chi - 1) / 2, (1 - chi^2), chi * (chi + 1) / 2];
                A_chi = Ne * A_x;
                
                Kea = Kea + (Be' * E * A_chi * Be) * (Je * wi);
            end

        otherwise
            error('El número de nodos especificado en A_x debe ser 2 o 3.');
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