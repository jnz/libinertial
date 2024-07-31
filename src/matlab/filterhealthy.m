function [healthy] = filterhealthy(x, P, chatty)
%FILTERHEALTH Check if there are issues in the current filter health.

healthy = true;

if all(all(isfinite(P))) ~= true
    healthy = false;
    if chatty
        fprintf('Matrix P is not finite\n');
    end
end

if all(isfinite(x)) ~= true
    healthy = false;
    if chatty
        fprintf('Vector x is not finite\n');
    end
end

if issymmetric(P) == false
    healthy = false;
    if chatty
        fprintf('Matrix P is not symmetric\n');
    end
end

try
    chol(P);
catch
    healthy = false;
    if chatty
        fprintf('Matrix P is not positive definite\n');
    end
end

end

