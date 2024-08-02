function [] = testcasegen()

kalman_udu_testcase();
kalman_udu_robust_testcase();
decorr_testcase();

end


function [] = kalman_udu_testcase()

R = eye(3)*0.25;
z = [16.2688, 17.9169, 16.8706]';
H = [ 8 1 6 1 ; 3 5 7 2; 4 9 2 3];
x = [1;1;1;1];

dz = z - H*x;

P = eye(4)*0.04;

[xexp, Pexp, chi2exp] = kalman_vanilla(x, P, dz, R, H);

[Uexp, dexp] = udu(Pexp);
R
z
H
x
xexp
Uexp'
dexp

end

function [] = kalman_udu_robust_testcase()

A = [1 1.2; 1.2 1];

x = [10; -5];
P = A*diag([0.1^2, 10.0^2])*A';
H = [1 -0.5; 0.1 5.0];
z = H*x - [0; 100]; % add outlier

B = [1.0 0.75; 0.2 5.0];
R = B*diag([0.25^2 1.5^2])*B';
dz = z - H*x;

[U, d] = udu(P);
chi2_threshold = 3.8415;
[x_robust_exp,U_exp,d_exp] = kalman_udu_robust(z,R,H,x,U,d,chi2_threshold,false);
P
H
x
z
R
x_robust_exp
U_exp'
d_exp

end

function [] = decorr_testcase()

x = [15; -2.5; 0];
H = [1 -0.5 0.25; 0.1 5.0 -2];
z = H*x; % add outlier

B = [1.0 0.75; 0.2 5.0];
R = B*diag([0.25^2 1.5^2])*B';

R
z
H

[G] = chol(R); % G'*G = R
zdecorr = (G')\z
Hdecorr = (G')\H

end

