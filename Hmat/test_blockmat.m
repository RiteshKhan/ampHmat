L=1;
d_dim = 2;
d_dimRootN = 80;
N = d_dimRootN^d_dim;
kernel_choice = 2;

depth = 3;
min_block_size = 100;
vareps = 1e-5;


% Particle location

m = ceil(N^(1/d_dim));
% 1D grid in each dimension
x = linspace(-L, L, m);
grids = cell(1, d_dim);
[grids{:}] = ndgrid(x);

% Total generated points
totalPts = m^d_dim;
pts = zeros(totalPts, d_dim);
for dim = 1:d_dim
    pts(:, dim) = grids{dim}(:);
end
% Trim down to exactly N points
particle_loc = pts(1:N, :);   % N x d

Ind = 1:N;
A = getMatrix(Ind, Ind, d_dim, particle_loc, kernel_choice);

subMat = A(1:N/256, N/256:N);

[U, S, V] = svd(subMat, 'econ');

rnk = sum(abs(diag(S)) > S(1,1) * vareps);

beta = norm(A);
err = norm(subMat-U(:,1:rnk)*S(1:rnk, 1:rnk)*V(:,1:rnk)')/beta;
fprintf("The norm 2 error: ");
disp(err);

fprintf("Simulation in different precisions\n\n");

sgm = S(1,1);
beta/sgm
prec = (vareps * beta)/sgm






function M = getMatrix(row_indices, col_indices, d_dim, particle_loc, kernel_choice)
    % row_indices, col_indices integer ROW vectors
    % Extract particle coordinates for the row and column sets
    r1 = particle_loc(row_indices, :);  % (n_rows x d_dim)
    r2 = particle_loc(col_indices, :);  % (n_cols x d_dim)

    % We'll build an (n_rows x n_cols) distance matrix
    n_rows = size(r1, 1);
    n_cols = size(r2, 1);
    M = zeros(n_rows, n_cols);

    % Compute pairwise distances (or kernel) in a fully vectorized manner
    %
    % For squared distances, do (r1 - r2)^2 dimension-wise, all at once:
    for dim = 1:d_dim
        diff = r1(:, dim) - r2(:, dim).';  % (n_rows x n_cols)
        M = M + diff.^2;
    end

    % Convert squared distance to actual distance
    M = sqrt(M);
    
    % (Matern kernel)
    if kernel_choice == 1
        M = exp(-M);
    end
    % 2D Laplacian -> log(R) with zero diagonal
    if kernel_choice == 2
        M = log(M);
        % Replace diagonal with 0
        idx = row_indices == col_indices.';
        M(idx) = 0;
    end
    
    % Example: 3D Laplacian -> 1/R with zero diagonal
    if kernel_choice == 3
        M = 1 ./ M;
        idx = row_indices == col_indices.';
        % Replace diagonal with 0
        M(idx) = 0;
    end
end