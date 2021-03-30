function Ez0 = gaussianSource(t, X, Y, E0, w0, t0, tp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%     R2 = X.^2+Y.^2;
%     Ez0 = E0 .* exp(-R2./w0^2) .* exp(-4*log(2)*(t-t0)^2/tp^2);
    xc = floor(length(X)/2);
    yc = floor(length(Y)/2);
    Ez0 = zeros(length(X),length(Y));
    Ez0(xc, yc) = E0*exp(-4*log(2)*(t-t0)^2/tp^2);
end

