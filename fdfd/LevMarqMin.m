function S = LevMarqMin(f, df, d2f, x, alpha, beta, mu, tol)
    n = 0;
    I = speye(length(x));
    dat = [n, f(x), norm(df(x))];

    while (norm(df(x)) > tol)
      dx = -(d2f(x) + mu*I) \ df(x);
      t = 1;
      while (f(x+t*dx) >= f(x) + alpha*(-dx')*(t*dx))
        t = beta * t;
      end
      x = x + t * dx;
      if isnan(x) == 1
          break;
      else
          n = n + 1;
          dat(n,:) = [n, f(x), norm(t*dx)];
      end
    end

    fprintf("%12s | %12s | %12s\n", ...
        "n", "f(x)", "||t*dx||")
    fprintf(repelem('-',40) + "\n")
    fprintf("%12d | %12g | %12g\n", ...
        dat')
    S.argmin = x;
    S.fmin = f(x);
    S.dnorm = norm(df(x));
end


