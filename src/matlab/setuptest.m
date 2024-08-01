function [] = setuptest()

rng(42);

epsilon = 1e-3;
epsilon_R = 1e-3;

H = [ 1e2 1 1; 1e2 1 1+epsilon];
fprintf('cond(H) = %.5e, rank(H) = %f\n', cond(H), rank(H));

P = eye(3);
x = [1e3; 0; 1e-03 ];

A = [1 0.5; 0 1];
% A = eye(2);
sigma = epsilon_R^2;
R = A*(eye(2)*sigma^2)*A';

S = H*P*H' + R
fprintf('cond(S) = %.5e, rank(S) = %f\n', cond(S), rank(S));

ztrue = H*x;
z = ztrue + A*randn(2,1)*sigma;

dz = z - H*x;

% Vanilla Test Numerical Stability

[x_vanilla, P_vanilla ] = kalman_vanilla(x, P, dz, R, H);

% Test Takasu Numerical Stability

[x_takasu, P_takasu ] = kalman_takasu(x, P, dz, R, H);

% UDU Test Numerical Stability

[U, d] = udu(P);

[x_udu,U_udu,d_udu] = kalman_udu(z,R,H,x,U,d);
P_udu = U_udu*diag(d_udu)*U_udu';

% Carlson Test Numerical Stability

Cxx_carlson = chol(P)';
[x_carlson,Cxx_carlson] = kalman_carlson(z,R,H,x,Cxx_carlson);
P_carlson = Cxx_carlson*Cxx_carlson';

% Check results

comparevec(x_vanilla, x_takasu, 1e-7, "x_vanilla", "x_takasu");
comparevec(x_udu,     x_takasu, 1e-7, "x_udu",     "x_takasu");
comparevec(x_carlson, x_takasu, 1e-7, "x_carlson", "x_takasu");
comparevec(x_carlson, x_udu,    1e-7, "x_carlson", "x_udu");

checkfilter(x_vanilla, P_vanilla, "Vanilla");
checkfilter(x_takasu,  P_takasu,  "Takasu");
checkfilter(x_carlson, P_carlson, "Carlson");
checkfilter(x_udu,     P_udu,     "UDU");

robust_udu_test();

end

function [] = comparevec(v1, v2, threshold, name1, name2)

diff = max(abs(v1-v2));
if (diff > threshold)
    fprintf(2, '%s differs from %s: %.3e\n', name1, name2, diff);
end

end

function [] = checkfilter(x, P, name)
chatty = true;
[healthy] = filterhealthy(x, P, chatty);
if (healthy == false)
    fprintf(2, 'Filter %s is not healthy\n', name);
end
end

function [] = robust_udu_test()

x = zeros(1,1);
P = eye(1)*0.1^2;
z = [-100];
H = [1];
R = eye(1)*0.1^2;
dz = z - H*x;

[U, d] = udu(P);
[x_non_robust,U_1,d_1] = kalman_udu(z,R,H,x,U,d);
[x_udu_non_robust,U_t,d_t] = kalman_udu_robust(z,R,H,x,U,d,0,false);
% in this case udu robust must equal the non-robust version:
assert(abs(x_non_robust-x_udu_non_robust) < 1e-8);

chi2_threshold = 3.8415;
[x_robust,U_2,d_2] = kalman_udu_robust(z,R,H,x,U,d,chi2_threshold,false);

[x_robust_takasu,P_takasu] = kalman_robust(x,P,dz,R,H,chi2_threshold,false);

assert(abs(x_robust-x_robust_takasu) < 1e-8);
assert(filterhealthy(x_robust,U_2*diag(d_2)*U_2'));
assert(filterhealthy(x_robust_takasu,P_takasu));

end

