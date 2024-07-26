function [] = benchmark()

rng(42);
epochs = 30000;

n = 4;
x = zeros(n,1); % Pos X, Pos Y, Vel X, Vel Y
P = diag([0.05^2, 0.05^2, 1.0^2, 1.0^2]);
Qnoise = diag([0.02^2, 0.02^2, 0.01^2, 0.01^2]);

H = [1 0 0 0; 0 1 0 0];
dt = 0.1;
Phi = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1 ];

sigma_pos = 0.05;
R = [sigma_pos^2, 0; 0, sigma_pos^2];

real_x = [0; 0; 2; 2];

tic
for i=1:epochs

    z = real_x(1:2) + randn(2, 1)*sigma_pos;
    dz = z - H*x;

    [x, P, chi2] = kalman_takasu(x, P, dz, R, H);

    % Next step:
    real_x = Phi*real_x;

    x = Phi*x;
    P = Phi*P*Phi' + Qnoise;
end
t = toc;
us_per_epoch = (t / epochs)*1e6;
fprintf('Total run-time: %.3f s (%.1f µs per epoch)\n', t, us_per_epoch);
fprintf('Expected:  %6.3f %6.3f %6.3f %6.3f\n', real_x(1), real_x(2), real_x(3), real_x(4));
fprintf('Estimated: %6.3f %6.3f %6.3f %6.3f\n', x(1), x(2), x(3), x(4));
error_x = real_x - x;
fprintf('Error:     %6.3f %6.3f %6.3f %6.3f\n', error_x(1), error_x(2), error_x(3), error_x(4));
fprintf('chi2: %.3f\n', chi2);
assert(max(abs(error_x(1:2))) < 2*sigma_pos);

% Only valid for rng(42) and 30000 epochs:
assert(abs(-0.039 - error_x(1)) < 1e-3);
assert(abs( 0.029 - error_x(2)) < 1e-3);
assert(abs(-0.017 - error_x(3)) < 1e-3);
assert(abs( 0.002 - error_x(4)) < 1e-3);

end
