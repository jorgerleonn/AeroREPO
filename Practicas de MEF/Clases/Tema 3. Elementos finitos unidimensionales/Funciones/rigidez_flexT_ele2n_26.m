function [Ke] = rigidez_flexT_ele2n_26(Le, A_x, I_x, ks, E, G, GL)
% RIGIDEZ_FLEXT_ELE2N_26 Matriz de rigidez a flexión/cortante de Timoshenko (2 nodos)
% Viga lagrangiana de 2 nodos con corrección por flexibilidad residual a flexión (RBF)
% para mitigar el bloqueo por cortante (shear locking).
%
% Grados de libertad locales: [w1, theta1, w2, theta2]'
%
% Sintaxis:
%   [Ke] = rigidez_flexT_ele2n_26(Le, A_x, I_x, ks, E, G, GL)
%
% Entradas:
%   Le   - Longitud del elemento [m o mm]
%   A_x  - Vector con el área en los 2 nodos: [A1; A2]
%   I_x  - Vector con el momento de inercia en los 2 nodos: [I1; I2]
%   ks   - Factor de corrección a cortante (ej. 5/6 para sección rectangular)
%   E    - Módulo de elasticidad longitudinal [Pa o MPa]
%   G    - Módulo de elasticidad transversal / cortante [Pa o MPa]
%   GL   - Orden de la cuadratura de Gauss-Legendre (1 a 5)
%
% Salida:
%   Ke   - Matriz de rigidez elemental (4x4) en coordenadas locales

    A_x = A_x(:);
    I_x = I_x(:);

    if length(A_x) ~= 2 || length(I_x) ~= 2
        error('Los vectores A_x e I_x deben contener exactamente 2 nodos.');
    end

    % Jacobiano 1D constante: dx/dchi = Le / 2
    Je = Le / 2;
    iJe = 1 / Je;

    Ke = zeros(4);
    [chi_ip, w_ip] = cuadratura(GL);

    for ip = 1:length(chi_ip)
        chi = chi_ip(ip);
        wi  = w_ip(ip);

        % Funciones de forma lagrangianas lineales
        Ne = [(1 - chi)/2, (1 + chi)/2];

        % Interpolación isoparamétrica de propiedades geométricas
        I_chi = Ne * I_x;
        A_chi = Ne * A_x;

        % Operadores de deformación cinemática (dof: [w1, theta1, w2, theta2]')
        % Curvatura: kappa = d(theta)/dx
        BeTf = [0, -1/2, 0, 1/2] * iJe;

        % Deformación a cortante: gamma = dw/dx - theta
        BeTs = [-1/2 * iJe, -(1 - chi)/2, 1/2 * iJe, -(1 + chi)/2];

        % Corrección por flexibilidad residual a cortante (Residual Bending Flexibility - RBF)
        % ks_mod = 1 / (1/ks + (A_chi * G * Le^2) / (12 * E * I_chi))
        ks_mod = 1 / ( (1 / ks) + (A_chi * G * Le^2) / (12 * E * I_chi) );

        % Integración numérica de rigidez de flexión + rigidez de cortante
        K_bending = (BeTf' * E * I_chi * BeTf) * (Je * wi);
        K_shear   = (BeTs' * (G * A_chi * ks_mod) * BeTs) * (Je * wi);

        Ke = Ke + K_bending + K_shear;
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