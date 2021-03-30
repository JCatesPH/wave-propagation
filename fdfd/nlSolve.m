function Xs = nlSolve(ez0, Ae, b, PNL, dPNL)
%nlSolve: Solves nonlinear system of Maxwell equations
%   Testing Newton method now
    F = @(x) sparse(Ae*x + PNL(x) - b);
	J = @(x) Ae + dPNL(x);
    
    res = 1;
    tol = 1e-9;
    ncalls = 0;
    while res > tol      
        ez1 = ez0 - J(ez0) \ F(ez0);
        res = norm(F(ez1));
        ez0 = ez1;
        ncalls = ncalls + 1;
    end
    Xs.ncalls = ncalls;
    Xs.sol = ez1;
    Xs.err = res;
end

