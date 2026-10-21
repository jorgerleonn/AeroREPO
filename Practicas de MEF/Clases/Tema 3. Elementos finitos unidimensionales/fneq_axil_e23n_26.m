function [Fneq_a] = fneq_axil_e23n_26(f_x, Le, GL)
% FNEQ_AXIL_E23N_26 Vector de fuerzas nodales equivalentes para carga axil distribuida
% Calcula el vector de fuerzas consistentes integrando f_eq = int(N' * q(x) * dx)
% en elementos barra 1D de 2 nodos (lineal) o 3 nodos (cuadrático).
%
% Sintaxis:
%   [Fneq_a] = fneq_axil_e23n_26(f_x, Le, GL)
%
% Entradas:
%   f_x  - Vector columna o fila con la carga distribuida en los nodos [N/m o N/mm]:
%          2 nodos: [q1; q2]
%          3 nodos: [q1; q2; q3] (nodo 2 es el punto medio)
%   Le   - Longitud del elemento [m o mm]
%   GL   - Orden de la cuadratura de Gauss-Legendre (1 a 4)
%
% Salida:
%   Fneq_a - Vector columna de fuerzas nodales equivalentes (2x1 o 3x1)

    f_x = f_x(:); % Asegurar vector columna
    num_nodos = length(f_x);

    % Jacobiano 1D constante: dx/dchi = Le / 2
    Je = Le / 2;

    [chi_ip, w_ip] = cuadratura(GL);
    Fneq_a = zeros(num_nodos, 1);

    switch num_nodos
        case 2
            % Elemento lineal de 2 nodos: chi in [-1, 1]
            for i = 1:length(chi_ip)
                chi = chi_ip(i);
                wi  = w_ip(i);

                Ne = [(1 - chi) / 2, (1 + chi) / 2];
                fx_chi = Ne * f_x; % Interpolación isoparamétrica de la carga

                Fneq_a = Fneq_a + (Ne' * fx_chi) * (Je * wi);
            end

        case 3
            % Elemento cuadrático de 3 nodos: chi in [-1, 1] (nodo medio en chi = 0)
            for i = 1:length(chi_ip)
                chi = chi_ip(i);
                wi  = w_ip(i);

                Ne = [chi * (chi - 1) / 2, (1 - chi^2), chi * (chi + 1) / 2];
                fx_chi = Ne * f_x; % Interpolación isoparamétrica de la carga

                Fneq_a = Fneq_a + (Ne' * fx_chi) * (Je * wi);
            end

        otherwise
            error('El vector f_x debe contener 2 o 3 valores nodales.');
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

        otherwise
            error('Orden de cuadratura no soportado. Usa un valor entre 1 y 4.');
    end
end