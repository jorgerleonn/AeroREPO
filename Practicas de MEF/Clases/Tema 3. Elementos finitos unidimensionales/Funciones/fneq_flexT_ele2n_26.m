function [Fneq_f] = fneq_flexT_ele2n_26(q_x, Le, GL)
% FNEQ_FLEXT_ELE2N_26 Vector de fuerzas nodales consistentes para viga de Timoshenko (2 nodos)
% Formulación lagrangiana lineal $C^0$ con grados de libertad locales:
%   d_e = [w1, theta1, w2, theta2]'
%
% Sintaxis:
%   [Fneq_f] = fneq_flexT_ele2n_26(q_x, Le, GL)
%
% Entradas:
%   q_x  - Vector con los valores de la carga transversal distribuida en los 2 nodos: [q1; q2] [N/m o N/mm]
%   Le   - Longitud del elemento [m o mm]
%   GL   - Orden de la cuadratura de Gauss-Legendre (1 a 4)
%
% Salida:
%   Fneq_f - Vector columna (4x1) de fuerzas y momentos nodales consistentes: [V1; M1; V2; M2]

    q_x = q_x(:); % Asegurar vector columna
    if length(q_x) ~= 2
        error('El vector de cargas q_x debe contener exactamente 2 nodos: [q1; q2].');
    end

    % Jacobiano 1D constante: dx/dchi = Le / 2
    Je = Le / 2;

    Fneq_f = zeros(4, 1);
    [chi_ip, w_ip] = cuadratura(GL);

    for ip = 1:length(chi_ip)
        chi = chi_ip(ip);
        wi  = w_ip(ip);

        % Funciones de forma lagrangianas lineales C^0
        N1 = (1 - chi) / 2;
        N2 = (1 + chi) / 2;

        % Interpolación de la carga transversal q(chi)
        qx_chi = [N1, N2] * q_x;

        % Matriz de interpolación respecto a [w1, theta1, w2, theta2]'
        % La carga transversal solo genera trabajo virtual sobre los desplazamientos w:
        % w(chi) = Nw * d_e -> Nw = [N1, 0, N2, 0]
        Nw = [N1; 0; N2; 0];

        % Integración consistente: int(Nw * q(x) * dx)
        Fneq_f = Fneq_f + (Nw * qx_chi) * (Je * wi);
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