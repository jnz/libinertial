Square Root Kalman Filter with Carlson Update and Schmidt-Householder Temporal Update
=====================================================================================

Jan Zwiener (jan@zwiener.org)

Start by a given state covariance matrix Qxx.
Perform

    Cxx = chol(Qxx)'; % note the transpose

where Cxx*Cxx' = Qxx;

Call

    [x,Cxx] = kalman_carlson(z,R,H,x,Cxx);

to update with a scalar measurement 'z'.
If you have multiple measurements, process them individually.
If the measurements are correlated, they are decorrelated inside
kalman_carlson(...) first with:

    [G] = chol(R); % G'*G = R
	zdecorr = (G')\z;
    Hdecorr = (G')\H;
    Rdecorr = eye(length(z));

The prediction step/temporal update of the Cxx matrix is done with the

    [Cxx] = kalman_carlson_predict(Cxx, phi, G, Cq)

Normal form:
    Qxx^(-) = phi*Qxx^(+)*phi' + G*Qnoise*G'

    Cq*Cq' = Qnoise -> from chol(...)' (note the transpose)


Example:

    % Init
    dt_sec = 1.0;
    x = [0 0 1]'; % state vector: position, vel, acceleration
    Qxx = diag([0.5 0.1 0.00001].^2); % initial state covariance
    Cxx = chol(Qxx)'; % Cxx*Cxx' = Qxx

    % Prediction
    phi = [1 dt_sec 0.5*dt_sec^2; 0 1 dt_sec; 0 0 1];
    x = phi*x; % state prediction
    % Add some noise to acceleration
    G = [0 0 dt_sec]';
    Cq = 0.001; % as std. dev. in m/s^2
    [Cxx] = kalman_carlson_predict(Cxx, phi, G, Cq); % predict Cxx

    % Fusion
    A = [1 0 0];
    z = 0.5; % position measurement
    R = 0.1^2; % measurement covariance
    [x,Cxx] = kalman_carlson(z,R,A,x,Cxx);

    % Calculate Qxx if required
    Qxx = Cxx*Cxx';
