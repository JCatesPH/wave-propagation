function Xs = nlSolve(ez0, Ae, b, PNL, dPNL, ddPNL, tol)
%nlSolve: Solves nonlinear system of Maxwell equations
%   Testing Newton method now
    F = @(x) sparse(Ae*x + PNL(x) - b); % Nonlinear function to solve
	J = @(x) Ae + dPNL(x); % Jacobian of the function
    H = @(x) ddPNL(x);
    
    conNL = normest(J(ez0)) / normest(F(ez0));
    fprintf("Nonlinear problem condition number: %e\n", conNL);
    fprintf('  n  |   con   |  res  |\n');
    fprintf('----------------------');
    
    % -- Lev-Marq Parameters --
    alpha = 0.1;
    beta = 0.5;
    
    res = 1;
    ncalls = 0;
    while res > tol    
        dez = -H(ez0) \ J(ez0);
        t = 1;
        while (F(ez0+t*dez) >= F(ez0) + alpha*(-dez')*(t*dez))
            t = beta * t;
        end
        ez1 = ez0 + t * dez;
        if isnan(ez1) == 1
            warning('Warning: NaN encountered in nonlinear solver.');
            break;
        else
            res = norm(F(ez1));
            ez0 = ez1;
            ncalls = ncalls + 1;
            conNL = normest(J(ez1)) / normest(F(ez1));
            fprintf(' %3d | %7.2e | %7.2e |\n', ncalls, conNL, res);
        end
    end
    Xs.ncalls = ncalls;
    Xs.sol = ez1;
    Xs.err = res;
end

