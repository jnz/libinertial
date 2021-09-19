function [x, P] = kalman_vanilla(x, P, dl, R, H)
%KALMAN_VANILLA Simple Kalman Filter implementation.
%
% x nx1 A priori state vector (size=n) at epoch k
% P nxn Covariance matrix of state vector x at epoch k
% dl mx1 Measurement difference vector (size=m) at epoch k
%        dl = measurement - predicted measurement
%        dl = y - H*x
% R mxm Covariance matrix of measurement vector y
% H mxn Observation matrix so that y = H*x
%
% Return value:
% x nx1 A posteriori state vector at epoch k (corrected by measurements)
% P nxn A posteriori covariance of x at epoch k
assert( ( (size(P,1)==size(P,2) ) && ... % symmetric covariance matrix
          (size(P,1)==size(x,1) ) && ...
          (size(x,2)==1 ) && ... % make sure x is a column vector
          (size(dl,2)==1 ) && ... % make sure dl is a column vector
          (size(dl,1)==size(R,1) ) && ...
          (size(R,1)==size(R,2) ) && ...
          (size(H,1)==size(dl,1) ) && ...
          (size(H,2)==size(x,1) ) && ...
          (nargin == 5) ), 'Invalid arguments');

% Add limits for the generated C-code. Otherwise there is no reason to limit
% the size here.
n_max = 32; % @satisfy{@req{3}}
m_max = 9;  % @satisfy{@req{3}}
assert(size(dl,1) <= m_max); 
assert(size(x,1) <= n_max);

% Vanilla:
K = (P*H')/(H*P*H' + R); % Kalman Gain Matrix K
dx = K*dl; % Correction to state by measurement
x = x + dx; % Update state vector, x is now the a posteriori state vec.
P = (eye(length(x)) - K*H)*P; % a posteriori covariance matrix
P = 0.5*(P+P');

end
