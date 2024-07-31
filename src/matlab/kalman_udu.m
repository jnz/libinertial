function [x,U,d] = kalman_udu(z,R,H,x,U,d)
% UDU' Bierman Filter Observation Step
%
% Inputs:
%  z   - measurement vector (m x 1)
%  R   - variance of measurement error (m x m)
%  H   - measurement sensitivity matrix (m x n)
%  x   - a priori estimate of state vector (n x 1)
%  U   - unit upper triangular factor of covariance matrix of a priori state uncertainty (n x n)
%  d   - diagonal vector with factor of covariance matrix of a priori state uncertainty (n x 1)
%
% Outputs:
%  x   - a posteriori estimate of state vector (n x 1)
%  U   - upper unit triangular UD factor of a posteriori state uncertainty covariance (n x n)
%  d   - diagonal UD factor vector of a posteriori state uncertainty covariance (n x 1)
%
% References:
%   1. Grewal, Weill, Andrews. "Global positioning systems, inertial
%      navigation, and integration". 1st ed. John Wiley & Sons, New York, 2001.

isdiagonal = isequal(R, diag(diag(R)));
assert(isdiagonal); % Measurements must not be correlated

% [G] = chol(R); % G'*G = R
% zdecorr = (G')\z;
% Hdecorr = (G')\H;

% Process measurements independently:
for i=1:size(H,1)
    [x,U,d] = kalman_udu_scalar(z(i),R(i,i),H(i,:),x,U,d);
end

end

function [x,U,d] = kalman_udu_scalar(z,R,H_line,xin,Uin,din)
% Inputs:
%  z   - scalar measurement (1 x 1)
%  R   - variance of measurement error (1 x 1)
%  H_line - row of measurement sensitivity matrix (m x 1)
%  x   - a priori estimate of state vector (n x 1)
%  U   - unit upper triangular factor of covariance matrix of a priori state uncertainty (n x n)
%  d   - diagonal vector with factor of covariance matrix of a priori state uncertainty (n x 1)
%
% Outputs:
%  x   - a posteriori estimate of state vector (n x 1)
%  U   - upper unit triangular UD factor of a posteriori state uncertainty covariance (n x n)
%  d   - diagonal UD factor vector of a posteriori state uncertainty covariance (n x 1)

x     = xin;      % Store inputs into outputs 
U     = Uin;      % because algorithm does in-place
d     = din;      % (destructive) calculation of outputs.
a     = U'*H_line';    % a is not modified, but
b     = din.*a;   % b is modified to become unscaled Kalman gain.
dz    = z - H_line*xin;
alpha = R;
gamma = 1/alpha;
for j=1:length(xin)
    beta   = alpha;
    alpha  = alpha + a(j)*b(j);
    lambda = -a(j)*gamma;
    gamma  = 1/alpha;
    d(j)   = beta*gamma*d(j);
    for i=1:j-1
        beta   = U(i,j);
        U(i,j) = beta + b(i)*lambda;
        b(i)   = b(i) + b(j)*beta;
    end
end
dzs = gamma*dz;  % apply scaling to innovations
x   = x + dzs*b; % multiply by unscaled Kalman gain

end

