#c = 3e8;
c = 1;
omeg = 1;
z = 1;

n1 = 1;
n2 = 1.5;
kz1 = sqrt(n1.^2 .* omeg.^2 ./c^2);
kz2 = sqrt(n2.^2 .* omeg.^2 ./c^2);

A1 = [1; 0]; 

M1 = [exp(-i*kz1*z), exp(i*kz1*z);
        -i*kz1.*exp(-i*kz1*z), i*kz1.*exp(-i*kz1*z)]
        
M2 = [exp(-i*kz2*z), exp(i*kz2*z);
        -i*kz2.*exp(-i*kz2*z), i*kz2.*exp(-i*kz2*z)]
        
A2 = M2 \ (M1 * A1)

S2 = abs(A2).^2

A1f = M1 \ (M2 * A2)