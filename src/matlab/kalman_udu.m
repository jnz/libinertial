function [x,U,d] = kalman_udu(z,R,H,x,U,d)
% UDU' Bierman Filter Observation Step
%
% Inputs:
%  z   - scalar measurement
%  R   - variance of measurement error
%  H   - measurement sensitivity (row) vector
%  z   - a priori estimate of state vector
%  U   - unit upper triangular factor of covariance matrix of a priori state uncertainty
%  d   - diagonal vector with factor of covariance matrix of a priori state uncertainty
%
% Outputs:
%  x   - a posteriori estimate of state vector
%  U   - upper triangular UD factor of a posteriori state uncertainty covariance
%  d   - diagonal UD factor vector of a posteriori state uncertainty covariance
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

function [x,U,d] = kalman_udu_scalar(z,R,H,xin,Uin,din)
x     = xin;      % Store inputs into outputs 
U     = Uin;      % because algorithm does in-place
d     = din;      % (destructive) calculation of outputs.
a     = U'*H';    % a is not modified, but
b     = din.*a;   % b is modified to become unscaled Kalman gain.
dz    = z - H*xin;
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

