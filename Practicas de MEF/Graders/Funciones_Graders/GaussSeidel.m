function [x, iterations] = GaussSeidel(A, b, tol, itmax)
    n = length(b);
    x = zeros(n, 1);
    normVal = 1e6;
    
    % * _*Tolerence for method*_
    if nargin == 2
        tol = 1e-8; 
        itmax = 2000;
    end
    
    iter = 0;
    while normVal > tol && iter < itmax
        x_old = x;
        for i = 1:n
            sigma = 0;
            for j = 1:i-1
                sigma = sigma + A(i,j)*x(j);
            end
            for j = i+1:n
                sigma = sigma + A(i,j)*x_old(j);
            end
            x(i) = (1/A(i,i))*(b(i) - sigma);
        end
        iter = iter + 1;
        normVal = norm(x_old - x);
    end
    iterations = iter;
end