function [x,Cxx] = kalman_carlson(z,R,H,x,Cxx)
% Carlson Cholesky filter observation step for linear systems
%
% Inputs:
%  z     - scalar measurement
%  R     - variance of measurement error
%  H     - measurement sensitivity matrix
%  x     - a priori estimate of state vector
%  Cxx   - upper triangular Cholesky factor of covariance matrix of a priori
%          state uncertainty
%
% Outputs:
%  x   - a posteriori estimate of state vector
%  Cxx - upper triangular Cholesky factor of covariance matrix of a posteriori
%        state uncertainty covariance
%
% References:
%   1. Grewal, Weill, Andrews. "Global positioning systems, inertial
%      navigation, and integration". 1st ed. John Wiley & Sons, New York, 2001.

isdiagonal = isequal(R, diag(diag(R)));
if (isdiagonal == false)
    [G] = chol(R); % G'*G = R
    zdecorr = (G')\z;
    Hdecorr = (G')\H;
    Rdecorr = eye(length(z));
else
    zdecorr = z;
    Hdecorr = H;
    Rdecorr = R;
end

% Process measurements independently:
for i=1:size(H,1)
    [x,Cxx] = kalman_carlson_scalar(zdecorr(i),Rdecorr(i,i),Hdecorr(i,:),x,Cxx);
end

end

function [x,Cxx] = kalman_carlson_scalar(z,R,H,xin,CxxIn)
C     = CxxIn;
alpha = R;
delta = z;
w = zeros(length(xin),1); % buffer

for j=1:length(xin)
    delta  = delta - H(j)*xin(j);
    sigma  = 0;
    for i=1:j
        sigma  = sigma + C(i,j)*H(i);
    end
    beta   = alpha;
    alpha  = alpha + sigma^2;
    gamma  = sqrt(alpha*beta);
    eta    = beta/gamma;
    zeta   = sigma/gamma;
    w(j)   = 0;
    for i=1:j
        tau    = C(i,j);
        C(i,j) = eta*C(i,j) - zeta*w(i);
        w(i)   = w(i) + tau*sigma;
    end
end
Cxx     = C;
epsilon = delta/alpha;
x       = xin + epsilon*w; % multiply by unscaled Kalman gain

end
