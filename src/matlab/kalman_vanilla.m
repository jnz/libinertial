function [x, P] = kalman_vanilla(x, P, dz, R, H)
%KALMAN_VANILLA Simple Kalman Filter implementation.
%
% x (n x 1) A priori state vector (size=n) at epoch k
% P (n x n) Covariance matrix of state vector x at epoch k
% dz (m x 1) Measurement difference vector (size=m) at epoch k
%        dz = measurement - predicted measurement
%        dz = z - H*x
% R (m x m) Covariance matrix of measurement vector z
% H (m x n) Observation matrix so that z = H*x
%
% Return value:
% x (n x 1) A posteriori state vector at epoch k (corrected by measurements)
% P (n x n) A posteriori covariance of x at epoch k
assert( ( (size(P,1)==size(P,2) ) && ... % symmetric covariance matrix
          (size(P,1)==size(x,1) ) && ...
          (size(x,2)==1 ) && ... % make sure x is a column vector
          (size(dz,2)==1 ) && ... % make sure dz is a column vector
          (size(dz,1)==size(R,1) ) && ...
          (size(R,1)==size(R,2) ) && ...
          (size(H,1)==size(dz,1) ) && ...
          (size(H,2)==size(x,1) ) && ...
          (nargin == 5) ), 'Invalid arguments');

K = (P*H')/(H*P*H' + R); % Kalman Gain Matrix K
dx = K*dz; % Correction to state by measurement
x = x + dx; % Update state vector, x is now the a posteriori state vec.
P = (eye(length(x)) - K*H)*P; % a posteriori covariance matrix
P = 0.5*(P+P');

end
