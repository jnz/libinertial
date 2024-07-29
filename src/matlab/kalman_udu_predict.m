function [x,U,D] = thornton(x,Phi,Uin,din,Gin,Q)
% UDU' Bierman-Thornton Filter Temporal / Prediction Step

%  Catherine Thornton's modified weighted Gram-Schmidt
%  orthogonalization method for the predictor update of
%  the U-D factors of the covariance matrix
%  of estimation uncertainty in Kalman filtering, see [1].
%
%  P^{+} = U*D*U' = Uin * diag(din) * Uin'
%
% Inputs:
%      x      - a posteriori state vector with size (n x 1)
%      Phi    - state transition matrix (n x n)
%      Uin    - unit upper triangular factor (U) of the modified Cholesky
%               factors (U-D factors) of the covariance matrix of
%               corrected state estimation uncertainty P^{+} (n x n)
%      din    - diagonal factor (d) vector (n x 1) of the U-D factors
%               of the covariance matrix of corrected estimation
%               uncertainty P^{+}, so that diag(din) = D.
%      Gin    - process noise distribution matrix (modified, if necessary to
%               make the associated process noise covariance diagonal) (n x r)
%      Q      - diagonal covariance matrix of process noise
%               in the stochastic system model (r x r)
% Outputs:
%      x      - predicted state vector (a priori)
%      U      - unit upper triangular factor (U) of the modified Cholesky
%               factors (U-D factors) of the covariance matrix of
%               predicted state estimation uncertainty P^{-}, so that
%               P^{-} = U*diag(d)*U'
%      d      - diagonal factor vector (d) of the U-D factors of the covariance
%               matrix of predicted estimation uncertainty (P-), so that
%               P^{-} = U*diag(d)*U'
%
% References:
%   1. Grewal, Weill, Andrews. "Global positioning systems, inertial
%      navigation, and integration". 1st ed. John Wiley & Sons, New York, 2001.

x     = Phi*x;
[n,r] = size(Gin); % get dimensions of state(n) and process noise (r)
G     = Gin;       % move to internal array for destructive updates
U     = eye(n);    % initialize lower triangular part of U
PhiU  = Phi*Uin;   % rows of [PhiU,G] are to be orthogonalized
for i=n:-1:1
    sigma = 0;
    for j=1:n
        sigma = sigma + PhiU(i,j)^2 *din(j);
        if (j <= r)
            sigma = sigma + G(i,j)^2 *Q(j,j);
        end;
    end;
    D(i,i) = sigma
    for j=1:i-1
        sigma = 0;
        for k=1:n
            sigma = sigma + PhiU(i,k)*din(k)*PhiU(j,k);
        end;
        for k=1:r
            sigma = sigma + G(i,k)*Q(k,k)*G(j,k);
        end
        U(j,i) = sigma/D(i,i);
        for k=1:n
            PhiU(j,k) = PhiU(j,k) - U(j,i)*PhiU(i,k);
        end
        for k=1:r
            G(j,k) = G(j,k) - U(j,i)*G(i,k);
        end
    end
end

