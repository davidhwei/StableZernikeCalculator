% A numerically stable routine for evaluating Zernike polynomials using
% the Janssen-Dirksen discrete Fourier transform method;
% Formulas:
%   U_n(r\cos\phi) = \sum_{m>=0}[1+(m>0)]R^m_n(r)\cos(m\phi);
%   R^m_n(r) = (1/2\pi)\int_{k=0}^{2\pi}U_n(r\cos\phi)\cos(m\phi)d\phi;
%   R^m_n(r) = (1/N)\sum_{k=0}^{N-1}U_n[r\cos(2\pi k/N)]\cos(2\pi mk/N);
%     U_n is the Chebyshev polynomial of the second kind;
%     U_n(x) = sin[(n+1)\arccos(x)]/sin[\arccos(x)];
%     0 <= r <= 1; m >= 0; n = m + 2p, p >= 0; N > m+n;
% Reference:
%   A. J.E.M. Janssen and P. Dirksen, "Computing Zernike polynomials of
% arbitrary degree using the discrete Fourier transform," J. Euro. Opt.
% Soc. - Rapid Pub. 2, 07012 (2007).
%
% [out] = R^m_n(r), with
%   m = AzimuthalOrder, scalar;
%   n = RadialOrder, scalar;
%   0 <= r <= 1, scalar;
%   [out] is a scalar when AzimuthalOrder is specified;
%   [out] is a vector when AzimuthalOrder is omitted;
%
% Copyright: David H. Wei, 2022
%
% Evaluate Zernike radial polynomials R^m_n(rho) for all m;
% where m and n are non-negative integers,
% m <= n and n-m is even;
% rho is in [0, 1], and can be a vector;

function [out] = ZernikeRn(n, rho)

if (size(rho, 1) > size(rho, 2))
    rVec = rho';
else
    rVec = rho;
end
Nr = length(rVec);

if (n < 1)
    out = ones(1, Nr);
else
    Nfft = 2^ceil(log2(2*n+1));
    kVec = cos((2*pi/Nfft)*(0:(Nfft-1)))';
    UnVec = ChebyshevU(n, kVec * rVec);
    ZernVec = real(ifft(UnVec, [], 1));
    Nm = 1 + floor(n/2);
    out = zeros(Nm, Nr);
    for im = 1:Nm
        m = 1 + mod(n,2) + 2*(im-1);
        out(im, :) = ZernVec(m, :);
    end
end
