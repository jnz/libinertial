function [x, P] = kalman_decorr(x, P, y, R, H)
%KALMAN_DECORR Kalman filter update routine based on measurement
%   decorrelation. Source section 8.1.3.1:
%   Grewal, Mohinder S., Lawrence R. Weill, and Angus P. Andrews. Global
%   positioning systems, inertial navigation, and integration. John Wiley &
%   Sons, 2007.
%
% x nx1 A priori state vector (size=n) at epoch k
% P nxn Covariance matrix of state vector x at epoch k
% y mx1 Measurement vector (size=m) at epoch k
% R mxm Covariance matrix of measurement vector y
% H mxn Observation matrix so that y = H*x
%
% Return value:
% x nx1 A posteriori state vector at epoch k (corrected by measurements)
% P nxn A posteriori covariance of x at epoch k
assert( ( (size(P,1)==size(P,2) ) && ... % symmetric covariance matrix
          (size(P,1)==size(x,1) ) && ...
          (size(x,2)==1 ) && ... % make sure x is a column vector
          (size(y,2)==1 ) && ... % make sure y is a column vector
          (size(y,1)==size(R,1) ) && ...
          (size(R,1)==size(R,2) ) && ...
          (size(H,1)==size(y,1) ) && ...
          (size(H,2)==size(x,1) ) && ...
          (nargin == 5) ), 'Invalid arguments');

% Add limits for the generated C-code. Otherwise there is no reason to limit
% the size here.
n_max = 32; % @satisfy{@req{3}}
m_max = 9;  % @satisfy{@req{3}}
assert(size(y,1) <= m_max);
assert(size(x,1) <= n_max);

[G] = chol(R); % G'*G = R

ydecorr = (G')\y;
Hdecorr = (G')\H;

for i=1:length(ydecorr)
    Hline = Hdecorr(i, :);

    % Vanilla form:
    % Rf = 1; % std. dev. is now exactly 1
    % K = P*Hline'/(Hline*P*Hline' + Rf);
    % dx = K*(ydecorr(i) - Hline*x);
    % x = x + dx;
    % P = P - K*Hline*P;

    % Grewal (2007): sec. 8.1.4 "Joseph stabilized implementation".
    K = (Hline*P*Hline' + 1)\(P*Hline');
    dx = K*(ydecorr(i) - Hline*x);
    x = x + dx;
    P = (eye(length(x)) - K*Hline)*P*(eye(length(x)) - K*Hline)' + K*K';
end

P = 0.5*(P + P');

end

