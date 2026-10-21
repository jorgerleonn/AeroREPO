function [Ke] = Ke_viga2D(E, Le, A, I)
    % rigideces escalares
    % E: módulo de Young
    % A: área de sección recta del elemento - constante
    % I: momento de inercia de la sección recta del elemento plano - constante
    % distancias relativas entre nodos locales 1 y 2 en el plano 2D
    EA = E*A; 
    EI = E*I;
    
    k11 = [ EA/Le,          0,          0;
                0, 12*EI/Le^3,  6*EI/Le^2; 
                0,  6*EI/Le^2,    4*EI/Le];
    k12 = [-EA/Le,          0,          0;
                0,-12*EI/Le^3,  6*EI/Le^2; 
                0, -6*EI/Le^2,    2*EI/Le];
    k22 = [ EA/Le,          0,          0;
                0, 12*EI/Le^3, -6*EI/Le^2; 
                0, -6*EI/Le^2,    4*EI/Le];
                
    k21 = k12';
    Ke = [k11 k12; k21 k22];
end
