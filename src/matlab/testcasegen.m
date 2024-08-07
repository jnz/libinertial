function [] = testcasegen()
% TESTCASEGEN This function helps to generate test cases and numbers for
% the C/C++ unit tests.

kalman_test();
kalman_udu_testcase();
kalman_udu_robust_testcase();
decorr_testcase();
kalman_takasu_robust();
thornton_test();

end

function [] = kalman_test()

fprintf('<default test>\n');
R = eye(3)*0.25
dz = [0.2688; 0.9169; -1.1294]
H = [8 1 6 1 ; 3 5 7 2; 4 9 2 3 ]
x = ones(4,1)
P = eye(4) * 0.04
[xnew, Pnew, chi2_takasu] = kalman_takasu(x, P, dz, R, H);

xnew
Pnew
fprintf('</default test>\n');

end

function [] = kalman_udu_testcase()

fprintf('<udu test>\n');
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

fprintf('</udu test>\n');

end

function [] = kalman_udu_robust_testcase()

fprintf('<robust udu test>\n');

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

fprintf('</robust udu test>\n');

end

function [] = decorr_testcase()

fprintf('<decorrelation test>\n');

x = [15; -2.5; 0];
H = [1 -0.5 0.25; 0.1 5.0 -2];
z = H*x;

B = [1.0 0.75; 0.2 5.0];
R = B*diag([0.25^2 1.5^2])*B';

R
z
H

[G] = chol(R); % G'*G = R
zdecorr = (G')\z
Hdecorr = (G')\H

fprintf('</decorrelation test>\n');

end

function [] = kalman_takasu_robust()

fprintf('<kalman_takasu outlier test>\n');

R = eye(3) * 0.25;
dz = [0.2688, 0.9169, -100.1294 ]';
H = [8 1 6 1 ; 3 5 7 2; 4 9 2 3 ];
x = ones(4,1);
P = eye(4) * 0.04;
[x, P, chi2_takasu] = kalman_takasu(x, P, dz, R, H);
chi2_takasu

fprintf('</kalman_takasu outlier test>\n');

end


function [] = thornton_test()

    Q = diag([0.1 0.2]);

    G = [1 0; 0 1; 0.5 0.5];

    B = [1 0.1 -0.2; 0.1 1.1 0.2; 0 0.2 1.0];
    P = B*eye(3)*B';

    x = [1; 2; 3];

    Phi = [1 0.5 0.25; 0 1 0.1; 0 0 1];

    [U,d] = udu(P);

    [x_exp,U_exp,d_exp] = kalman_udu_predict(x,Phi,U,d,G,Q);

    x_ref = Phi*x;
    P_ref = Phi*P*Phi' + G*Q*G';
    [U_ref,d_ref] = udu(P_ref);
    assert(max(max(abs(U_ref-U_exp))) < 1e-10);

    k=0;

end

