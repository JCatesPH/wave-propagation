function S = BFGS(f, df, x, Ck, alpha, beta, tol)
    n = 0;
    I = speye(length(x));  
    dat = [n, f(x), norm(Ck*df(x))];

    while (norm(df(x)) > tol)
      df0 = df(x);
      dx = -Ck * df0;
      t = 1;
      while (f(x+t*dx) >= f(x) + alpha * df0' * (t*dx))
        t = beta * t;
        if (norm(t*dx) < eps)
            t = beta;
            break;
        end
      end
      dn = t * dx;
      x = x + dn;
      gn = df(x) - df0;
      dngnin = gn' * dn;
      Ck = (I-(dn*gn')/dngnin)* Ck *(I-(gn*dn')/dngnin) + (dn*dn')/dngnin;
      if isnan(x) == 1
          break;
      else
          n = n + 1;
          dat(n,:) = [n, f(x), norm(dn)];
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
    S.hist = dat;
    S.xhist = xmat;
end

