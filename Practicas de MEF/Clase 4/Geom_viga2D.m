function [Le, Te] = Geom_viga2D(n1, n2, DIR)
    % n1, n2: coordenadas globales de los nodos locales 1 y 2
    % base global de referencia
    n_X = DIR(:,1); 
    n_Y = DIR(:,2); 
    n_Z = DIR(:,3);
    
    % diferencias de longitud en coordenadas globales
    delta = n2 - n1; 
    delta_x = delta(1); 
    delta_y = delta(2); % distancia entre nodos
    Le = sqrt(delta_x^2 + delta_y^2); % Le longitud global del elemento
    
    C = delta_x/Le; % coseno del ángulo en el plano
    S = delta_y/Le; % seno del ángulo en el plano
    
    % vectores de la base local x'y' en coordenadas globales
    n_x = [C; S; 0]; 
    n_y = [-S; C; 0]; 
    n_z = [0; 0; 1];
    
    % Re matriz de cambio de base global-local
    CXx = n_x'*n_X; CYx = n_x'*n_Y; CZx = n_x'*n_Z;
    CXy = n_y'*n_X; CYy = n_y'*n_Y; CZy = n_y'*n_Z;
    CXz = n_z'*n_X; CYz = n_z'*n_Y; CZz = n_z'*n_Z;
    
    Re = [CXx CYx CZx; CXy CYy CZy; CXz CYz CZz];
    % Re = [C S 0; -S C 0; 0 0 1]; % escrito directamente sin operaciones
    
    R0 = zeros(3);
    Te = [Re R0; R0 Re]; % Te matriz de transformación global-local de elemento de 2 nodos
end