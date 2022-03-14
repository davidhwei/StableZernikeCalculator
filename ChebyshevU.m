% Evaluate Chebyshev polynomial of the second kind U_n(x);
% where n is a non-negative integer, x is in [-1, 1];

function y = ChebyshevU(n, x)

format long

AlmostZero = 1e-7;

phi = acos(x); % phi is in [0, pi]
sin_phi_raw = sin(phi);
phi_big = abs(sin_phi_raw) > AlmostZero;
sin_phi = sin_phi_raw .* phi_big + AlmostZero * (1 - phi_big);
sin_np1_phi = sin((n+1)*phi) .* phi_big + ((n+1)*(1-2*(x<0)*mod(n,2))*AlmostZero) .* (1 - phi_big);
y = sin_np1_phi ./ sin_phi;
