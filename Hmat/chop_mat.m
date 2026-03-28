%% Chop low-rank factors
function [U, V] = chop_mat(U, S, V)
    global opt;
    U = mchop(U);
    V = diag(S) * mchop(V');
end