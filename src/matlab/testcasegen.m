function [] = testcasegen()

kalman_udu_testcase();

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

