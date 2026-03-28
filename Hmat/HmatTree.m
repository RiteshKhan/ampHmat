classdef HmatTree < handle
    properties
        K                  % Kernel object (handle or struct)
        nLevels            % Number of levels in the tree
        hodlr_levels       % Number of levels where HODLR is performed
        N                  % Number of particles
        N_max              % Max size of N for which we build the entire matrix

        kernel_choice      % Choice of kernel function
        is_sym             % Is kernel matrix symmetric?
        isNBDcompressed    % Note, it decides whether the nbd will be compressed at joining level

        L                  % Semi-length of the simulation box
        smallestBoxSize    % L / 2^nLevels

        nBoxesPerLevel     % 1×(nLevels+1) array of number of boxes per level
        boxRadius          % 1×(nLevels+1) array of box radii

        tree               % Cell array {level}{box} storing Hmat_node objects

        LR_epsilon         % Tolerance for low-rank compression
        nParticlesInLeafAlong1D
        nParticlesInLeaf

        Nodes1D            % 1D node coordinates
        Nodes              % Full list of particle positions (cell or matrix)
        gridPoints

        pow2d
        eta
        d_dim
        sqrt_d_dim

        particle_loc

        computed_potential % Column vector of size N
        collective_charge  % Column vector of size N

        prec_settings
        working_prec

        LR_storage
        Dense_storage
        storage

        normOrder
        sq_matrixNorm  % frob(A) * frob(A)

        Hmat_storage
        HODLR_strorage
        nadm_storage;

        sortIdx
        unitRoundOff {mustBeNonNan, mustBeFinite, mustBeNumeric}

        Hmat
    end

    methods
        function obj = HmatTree(d_dim, eta, N, N_max, nLevels, hodlr_levels, nParticlesInLeafAlong1D, L, LR_epsilon, kernel_choice, is_sym, isNBDcompressed, prec, u)
            obj.d_dim = d_dim;
            obj.sqrt_d_dim = sqrt(d_dim);
            obj.eta = eta;

            obj.pow2d = 2^d_dim;
            obj.nLevels = nLevels;
            obj.hodlr_levels = hodlr_levels;
            obj.L = L;
            obj.nParticlesInLeafAlong1D = nParticlesInLeafAlong1D;
            obj.nParticlesInLeaf = nParticlesInLeafAlong1D^(obj.pow2d);
            obj.LR_epsilon = LR_epsilon;
            obj.N = N;
            obj.N_max = N_max;

            obj.kernel_choice = kernel_choice;
            obj.is_sym = is_sym;
            obj.isNBDcompressed = isNBDcompressed;

            obj.working_prec = u;
            obj.prec_settings = [prec_chain(u), prec];
            [obj.unitRoundOff, obj.sortIdx] = sort_by_u(obj, obj.prec_settings);

            % Initialize box counts and radii
            obj.nBoxesPerLevel = zeros(1, nLevels+1);
            obj.boxRadius      = zeros(1, nLevels+1);

            obj.nBoxesPerLevel(1) = 1;
            obj.boxRadius(1) = L;
            for k = 2:nLevels+1
                obj.nBoxesPerLevel(k) = obj.pow2d * obj.nBoxesPerLevel(k-1);
                obj.boxRadius(k) = obj.boxRadius(k-1) / 2;
            end

            obj.smallestBoxSize = obj.boxRadius(end);

            % Initialize tree storage as a cell array of levels
            obj.tree = cell(1, nLevels+1);
        end


        function createTree(obj)
            % --- Create root box (level 1 in MATLAB, level 0 in C++) ---
            root = Hmat_node(obj.d_dim);
            root.boxNumber = 1;           % Parent of root level is 0
            root.parentNumber = 0;

            % Assign children numbers 0 to pow2d-1
            root.childrenNumbers = 1:obj.pow2d;

            % Store in tree{1}{1}
            obj.tree{1} = {root};  % Level 1 contains one box

            % --- Create remaining levels ---
            for j = 2:(obj.nLevels+1)  % MATLAB level index (C++ j = 1 to nLevels)
                nBoxes = obj.nBoxesPerLevel(j);
                levelBoxes = cell(1, nBoxes);
                for k = 1:nBoxes
                    box = Hmat_node(obj.d_dim);
                    box.boxNumber = k;
                    box.parentNumber = ceil(k / obj.pow2d);        % 1-based parent
                    box.childrenNumbers = obj.pow2d*(k-1) + (1:obj.pow2d);  % 1-based children
                    levelBoxes{k} = box;
                end
                obj.tree{j} = levelBoxes;
            end
        end


        function assign_Center_Locations(obj)
            % Set the root center once
            obj.tree{1}{1}.center = zeros(1, obj.d_dim);

            % Precompute binary child offset patterns (±1 per dimension)
            % Rows correspond to children, columns to dimensions
            binMat = dec2bin(0:obj.pow2d-1, obj.d_dim) == '1';
            signMat = double(binMat) * 2 - 1;  % 1 for bit=1, -1 for bit=0
            % Loop over levels
            for j = 1:obj.nLevels
                J = j + 1;
                shift = 0.5 * obj.boxRadius(j);

                nBoxes = obj.nBoxesPerLevel(j);
                for k = 1:nBoxes
                    parentCenter = obj.tree{j}{k}.center;

                    % Assign child centers
                    for c = 1:obj.pow2d
                        idxChild = obj.pow2d*(k-1) + c;
                        obj.tree{J}{idxChild}.center = parentCenter + shift * signMat(c, :);
                    end
                end
            end
        end


        % --- Euclidean distance ---
        function d = get_Euclid_dist(~, p1, p2)
            diff = p1 - p2;
            d = sqrt(sum(diff.^2));
        end
        

        % --- Admissibility condition (Euclidian norm) ---
        function tf = admCondition(obj, j, k, l)
            c1 = obj.tree{j}{k}.center;
            c2 = obj.tree{j}{l}.center;

            dist = obj.get_Euclid_dist(c1, c2);
            lhs = obj.sqrt_d_dim * 2 * obj.boxRadius(j);
            rhs = obj.eta * dist;

            tf = (lhs - rhs) <= 1e-16;
        end


        %  % --- Admissibility condition (Max norm) ---
        % function tf = admCondition(obj, j, k, l)
        %     c1 = obj.tree{j}{k}.center;
        %     c2 = obj.tree{j}{l}.center;
        %     tf = false;
        %     for dim = 1:obj.d_dim
        %         if 2 * obj.boxRadius(j) - obj.eta * abs(c1(dim)-c2(dim)) < 1e-16
        %             tf = true;
        %             return;
        %         end
        %     end
        % end


        function assign_Box_Interactions(obj, j, k)
            % ----- Parent level and box (1-based indexing) -----
            pj = j - 1;                  % parent level
            pk = ceil(k / obj.pow2d);       % parent box index at level pj

            % ----- 1. Children of parent's neighbors -----
            parentNeighbors = obj.tree{pj}{pk}.neighborNumbers;

            for idx = 1:numel(parentNeighbors)
                pn = parentNeighbors(idx);  % parent neighbor (1-based)

                for i = 1:obj.pow2d
                    cpn = (pn-1)*obj.pow2d + i;  % child of pn (1-based)

                    if obj.admCondition(j, k, cpn)
                        obj.tree{j}{k}.interactionList(end+1) = cpn;
                    else
                        obj.tree{j}{k}.neighborNumbers(end+1) = cpn;
                    end
                end
            end
            % ----- 2. Siblings -----
            for i = 1:obj.pow2d
                cpk = (pk-1)*obj.pow2d + i;  % sibling box index in same level (1-based)

                if k ~= cpk
                    if obj.admCondition(j, k, cpk)
                        obj.tree{j}{k}.interactionList(end+1) = cpk;
                    else
                        obj.tree{j}{k}.neighborNumbers(end+1) = cpk;
                    end
                end
            end
        end


        function assign_Tree_Interactions(obj)
            % Assign interactions for all levels except the root (level 1)
            for j = 2:(obj.nLevels-obj.hodlr_levels+1)
                % Assign interactions for all boxes at level j
                nBoxes = obj.nBoxesPerLevel(j);
                for k = 1:nBoxes
                    obj.assign_Box_Interactions(j, k);
                end
            end
        end


        % function printIL(obj)
        %     maxIL = 0;
        %     for j = 1:(obj.nLevels+1)  % MATLAB levels = 1..nLevels+1
        %         nBoxes = obj.nBoxesPerLevel(j);
        %         for k = 1:nBoxes
        %             fprintf('Level (j): %d   Box (k): %d\n', j, k);
        %
        %             ilist = obj.tree{j}{k}.interactionList;
        %             if ~isempty(ilist)
        %                 fprintf('%s\n', sprintf('%d,', ilist));
        %                 maxIL = max(maxIL, numel(ilist));
        %             else
        %                 fprintf('No interactions\n');
        %             end
        %             fprintf('------------------------\n');
        %         end
        %     end
        %     fprintf('Maximum interaction list size: %d\n', maxIL);
        % end


        % function printN(obj)
        %     maxN = 0;
        %     for j = 1:(obj.nLevels+1)
        %         nBoxes = obj.nBoxesPerLevel(j);
        %         for k = 1:nBoxes
        %             fprintf('Level (j): %d   Box (k): %d\n', j, k);
        %
        %             nlist = obj.tree{j}{k}.neighborNumbers;
        %             if ~isempty(nlist)
        %                 fprintf('%s\n', sprintf('%d,', nlist));
        %                 maxN = max(maxN, numel(nlist));
        %             else
        %                 fprintf('No neighbors\n');
        %             end
        %             fprintf('------------------------\n');
        %         end
        %     end
        %     fprintf('Maximum neighbor list size: %d\n', maxN);
        % end


        function assign_Particle_Locations(obj, choice)  % Uniform or Random
            if strcmp(choice, 'uniform')
                m = ceil(obj.N^(1/obj.d_dim));
                x = linspace(-obj.L, obj.L, m);    % uniform grid

                grids = cell(1, obj.d_dim);
                [grids{:}] = ndgrid(x);

                totalPts = m^obj.d_dim;
                pts = zeros(totalPts, obj.d_dim);
                for dim = 1:obj.d_dim
                    pts(:, dim) = grids{dim}(:);
                end

                obj.particle_loc = pts(1:obj.N, :);

            elseif strcmp(choice, 'random')
                % Uniform random points in [-L, L]^d
                obj.particle_loc = (2*rand(obj.N, obj.d_dim) - 1) * obj.L;
            else
                error("Please select proper particle location: 'uniform' or 'random'");
            end
        end


        function assign_Particle_Locations_Surface(obj, choice)
            % Generate N points uniformly on the surface of a sphere (or circle)
            % inside [-L, L]^d.  d = obj.d_dim ∈ {2,3}.

            d = obj.d_dim;

            if d == 2
                if strcmp(choice, 'circle')
                    theta = 2*pi*rand(obj.N,1);
                    %theta = sort(theta);
                    r = obj.L;            % radius = L
                    x = r * cos(theta);
                    y = r * sin(theta);
                    obj.particle_loc = [x y];
                elseif strcmp(choice, 'star')
                    n_tips = 3;               % number of star tips
                    theta = linspace(0, 2*pi, obj.N + 1)';
                    theta(end) = [];          % remove last point to avoid duplication
                    r = obj.L * (0.5 + 0.5 * sin(n_tips*theta).^2);  % radial modulation
                    x = r .* cos(theta);
                    y = r .* sin(theta);
                    obj.particle_loc = [x y];
                else
                    error("Please select proper particle location!!!");
                end

            elseif d == 3
                if strcmp(choice, 'sphere')

                    % Fibonacci lattice
                    i = (0:obj.N-1)';
                    phi = acos(1 - 2*(i+0.5)/obj.N);
                    theta = pi*(1 + sqrt(5))*(i+0.5);

                    r = obj.L;

                    x = r .* sin(phi) .* cos(theta);
                    y = r .* sin(phi) .* sin(theta);
                    z = r .* cos(phi);

                    obj.particle_loc = [x y z];
                elseif strcmp(choice, 'polyhedron')
                    theta = 2*pi*rand(obj.N,1);
                    phi   = acos(2*rand(obj.N,1)-1);
                    k_theta = 3;   % spikes around equator
                    k_phi   = 4;   % spikes from pole to pole
                    a       = 0.5; % spike amplitude
                    r = obj.L * (1 + a*cos(k_theta*theta).*cos(k_phi*phi)) / (1 + a);
                    x = r .* sin(phi) .* cos(theta);
                    y = r .* sin(phi) .* sin(theta);
                    z = r .* cos(phi);
                    obj.particle_loc = [x y z];
                else
                    error("Please select proper particle location!!!");
                end
            else
                error("assign_Particle_Locations_Sphere supports only d=2 or d=3.");
            end
        end


        function print_particle_location(obj)
            disp(obj.particle_loc);
        end


        function assign_Charge_Locations(obj)
            % Assign all particle indices to the root box initially
            obj.tree{1}{1}.chargeLocations = 1:obj.N;
            % Loop through each level to push particles down to children
            for j = 1:obj.nLevels
                for k = 1:obj.nBoxesPerLevel(j)
                    parentCL = obj.tree{j}{k}.chargeLocations;
                    J  = j + 1;                 % Parent level
                    Kp = (k - 1) * obj.pow2d + 1;   % Starting index for children at level J
                    for idx = 1:numel(parentCL)
                        pIndex = parentCL(idx);    % Particle index (1-based)
                        center = obj.tree{j}{k}.center;

                        % Determine binary code for child index
                        bits = false(1, obj.d_dim);
                        for dim = 1:obj.d_dim
                            if obj.particle_loc(pIndex, dim) <= center(dim)
                                bits(dim) = 0;
                            else
                                bits(dim) = 1;
                            end
                        end

                        c = bin2dec(char('0' + bits));        % Convert bits -> integer
                        childBoxIndex = Kp + c;

                        % Append to that child's charge locations
                        obj.tree{J}{childBoxIndex}.chargeLocations(end + 1 ) = pIndex;
                    end
                end
            end
        end


        function assign_NonLeaf_Charge_Locations(obj)
            % Move upward from level nLevels to 2
            for j = obj.nLevels:-1:2
                for k = 1:obj.nBoxesPerLevel(j)
                    % Clear current node's charge list
                    obj.tree{j}{k}.chargeLocations = [];

                    % For each child, append its chargeLocations
                    for c = 1:obj.pow2d
                        childIndex = (k - 1)*obj.pow2d + c;
                        obj.tree{j}{k}.chargeLocations = [obj.tree{j}{k}.chargeLocations,obj.tree{j+1}{childIndex}.chargeLocations];
                    end
                end
            end
        end


        function print_charge_locations(obj)
            for j = 1:(obj.nLevels + 1)
                % Assign interactions for all boxes at level j
                nBoxes = obj.nBoxesPerLevel(j);
                for k = 1:nBoxes
                    fprintf("Box %d \n", k);
                    disp(obj.tree{j}{k}.chargeLocations);
                end
                fprintf("\n======== ** ========\n");
            end
        end

        %% Assign per-node charges vector using global charges (Nx1)
        function assign_Charges(obj, vec)
            obj.collective_charge = vec(:); % store column vector
            % Loop levels 2..end (leaf level included)
            for j = 2:(obj.nLevels+1)
                for k = 1:obj.nBoxesPerLevel(j)
                    locs = obj.tree{j}{k}.chargeLocations;
                    if isempty(locs)
                        obj.tree{j}{k}.charges = [];
                    else
                        % locs are indices into charges (1-based)
                        obj.tree{j}{k}.charges = obj.collective_charge(locs);
                    end
                end
            end
        end


        function M = getMatrix(obj, row_indices, col_indices)
            % row_indices, col_indices integer ROW vectors
            % Extract particle coordinates for the row and column sets
            r1 = obj.particle_loc(row_indices, :);  % (n_rows x d_dim)
            r2 = obj.particle_loc(col_indices, :);  % (n_cols x d_dim)

            % We'll build an (n_rows x n_cols) distance matrix
            n_rows = size(r1, 1);
            n_cols = size(r2, 1);
            M = zeros(n_rows, n_cols);

            % Compute pairwise distances (or kernel) in a fully vectorized manner

            for dim = 1:obj.d_dim
                diff = r1(:, dim) - r2(:, dim).';  % (n_rows x n_cols)
                M = M + diff.^2;
            end

            switch obj.kernel_choice
                % Matern kernel
                case 1
                    M = exp(-sqrt(M));

                    % 2D Laplacian -> log(R) with zero diagonal
                case 2
                    M = 0.5*log(M); % log(sqrt(M))
                    % Replace diagonal with 0
                    idx = bsxfun(@eq, row_indices(:), col_indices(:).');
                    M(idx) = 0;

                    % 3D Laplacian -> 1/R with zero diagonal
                case 3
                    M = 1 ./ sqrt(M);
                    idx = bsxfun(@eq, row_indices(:), col_indices(:).');
                    % Replace diagonal with 0
                    M(idx) = 0;

                    % Gaussian kernel
                case 4
                    M = exp(-M/2);

                    % 1/r^2
                case 5
                    M = 1 ./ M;
                    idx = bsxfun(@eq, row_indices(:), col_indices(:).');
                    % Replace diagonal with 0
                    M(idx) = 0;

                otherwise
                    error('Unknown kernel choice');
            end
        end


        function print_Matrix(obj)
            M = obj.getMatrix(obj.tree{1}{1}.chargeLocations, obj.tree{1}{1}.chargeLocations);
            display(M);
        end


        function compress_Hmat_Block_Matrices(obj)
            obj.sq_matrixNorm = 0;
            obj.Hmat_storage = zeros(obj.nLevels-obj.hodlr_levels+1, 1);
            for j = 2:(obj.nLevels-obj.hodlr_levels+1)
                nBoxes = obj.nBoxesPerLevel(j);
                for k = 1:nBoxes
                    boxK = obj.tree{j}{k};
                    I_list = boxK.interactionList;
                    boxK.U = cell(numel(I_list),1);
                    boxK.S = cell(numel(I_list),1);
                    boxK.V = cell(numel(I_list),1);
                    boxK.sq_norm_block_mat = cell(numel(I_list),1);
                    for idx = 1:numel(I_list)
                        ki = I_list(idx);
                        boxKI = obj.tree{j}{ki};
                        if(obj.is_sym)
                            if(~boxK.isVisited(ki))
                                block_mat = obj.getMatrix(boxK.chargeLocations, boxKI.chargeLocations);
                                if isempty(block_mat)
                                    boxK.U{ki} = [];
                                    boxK.S{ki} = [];
                                    boxK.V{ki} = [];
                                    boxK.sq_norm_block_mat{ki} = 0;
                                    continue
                                end

                                % Economy-size SVD
                                [boxK.U{ki}, boxK.S{ki}, boxK.V{ki}, ~, boxK.sq_norm_block_mat{ki}] = mp_compress_m(block_mat, "rsvd", obj.LR_epsilon); % use the norm to get xi and set precs

                                obj.sq_matrixNorm = obj.sq_matrixNorm + 2 * boxK.sq_norm_block_mat{ki};
                                boxKI.isVisited(k) = true;
                            end
                        else
                            block_mat = obj.getMatrix(boxK.chargeLocations, boxKI.chargeLocations);
                            if isempty(block_mat)
                                boxK.U{ki} = [];
                                boxK.S{ki} = [];
                                boxK.V{ki} = [];
                                boxK.sq_norm_block_mat{ki} = 0;
                                continue
                            end
                            % Economy-size SVD
                            [boxK.U{ki}, boxK.S{ki}, boxK.V{ki}, ~, boxK.sq_norm_block_mat{ki}] = mp_compress_m(block_mat, "rsvd", obj.LR_epsilon); % use the norm to get xi and set precs
                            obj.sq_matrixNorm = obj.sq_matrixNorm + boxK.sq_norm_block_mat{ki};
                        end
                    end
                end
            end
        end


        function compress_Junction_Block_Matrices(obj)
            j = obj.nLevels-obj.hodlr_levels+1;
            nBoxes = obj.nBoxesPerLevel(j);
            obj.nadm_storage = 0;
            for k = 1:nBoxes
                boxK = obj.tree{j}{k};
                for idx = 1:numel(boxK.neighborNumbers)
                    nn = boxK.neighborNumbers(idx);
                    boxNN = obj.tree{j}{nn};
                    if(obj.is_sym)
                        if(~boxK.isVisited(nn))
                            block_mat = obj.getMatrix(boxK.chargeLocations, boxNN.chargeLocations);
                            if isempty(block_mat)
                                boxK.U{nn} = [];
                                boxK.S{nn} = [];
                                boxK.V{nn} = [];
                                boxK.sq_norm_block_mat{nn} = 0;
                                continue
                            end

                            % Economy-size SVD
                            [boxK.U{nn}, boxK.S{nn}, boxK.V{nn}, ~, boxK.sq_norm_block_mat{nn}] = mp_compress_m(block_mat, "rsvd", obj.LR_epsilon); % use the norm to get xi and set precs

                            obj.sq_matrixNorm = obj.sq_matrixNorm + 2 * boxK.sq_norm_block_mat{nn};
                            boxNN.isVisited(k) = true;
                        end
                    else
                        block_mat = obj.getMatrix(boxK.chargeLocations, boxNN.chargeLocations);
                        if isempty(block_mat)
                            boxK.U{nn} = [];
                            boxK.S{nn} = [];
                            boxK.V{nn} = [];
                            boxK.sq_norm_block_mat{nn} = 0;
                            continue
                        end
                        % Economy-size SVD
                        [boxK.U{nn}, boxK.S{nn}, boxK.V{nn}, ~, boxK.sq_norm_block_mat{nn}] = mp_compress_m(block_mat, "rsvd", obj.LR_epsilon); % use the norm to get xi and set precs
                        obj.sq_matrixNorm = obj.sq_matrixNorm + boxK.sq_norm_block_mat{nn};

                    end
                end
            end
        end


        function compress_HODLR_Block_Matrices(obj)
            obj.HODLR_strorage = zeros(obj.nLevels+1, 1);
            for j = (obj.nLevels-obj.hodlr_levels+2):(obj.nLevels+1)
                nBoxes = obj.nBoxesPerLevel(j);
                for k = 1:nBoxes
                    boxK = obj.tree{j}{k};
                    boxK.U = cell(numel(obj.pow2d),1);
                    boxK.S = cell(numel(obj.pow2d),1);
                    boxK.V = cell(numel(obj.pow2d),1);
                    boxK.sq_norm_block_mat = cell(numel(obj.pow2d),1);
                    pN = ceil(k / obj.pow2d);
                    for c = 1:obj.pow2d
                        ki = obj.tree{j-1}{pN}.childrenNumbers(c);
                        if(ki ~= k)
                            boxKI = obj.tree{j}{ki};
                            if(obj.is_sym)
                                if(~boxK.isVisited(ki))
                                    block_mat = obj.getMatrix(boxK.chargeLocations, boxKI.chargeLocations);
                                    if isempty(block_mat)
                                        boxK.U{ki} = [];
                                        boxK.S{ki} = [];
                                        boxK.V{ki} = [];
                                        boxK.sq_norm_block_mat{ki} = 0;
                                        continue
                                    end

                                    % Economy-size SVD
                                    [boxK.U{ki}, boxK.S{ki}, boxK.V{ki}, ~, boxK.sq_norm_block_mat{ki}] = mp_compress_m(block_mat, "rsvd", obj.LR_epsilon); % use the norm to get xi and set precs

                                    obj.sq_matrixNorm = obj.sq_matrixNorm + 2 * boxK.sq_norm_block_mat{ki};
                                    boxKI.isVisited(k) = true;
                                end
                            else
                                block_mat = obj.getMatrix(boxK.chargeLocations, boxKI.chargeLocations);
                                if isempty(block_mat)
                                    boxK.U{ki} = [];
                                    boxK.S{ki} = [];
                                    boxK.V{ki} = [];
                                    boxK.sq_norm_block_mat{ki} = 0;
                                    continue
                                end
                                % Economy-size SVD
                                [boxK.U{ki}, boxK.S{ki}, boxK.V{ki}, ~, boxK.sq_norm_block_mat{ki}] = mp_compress_m(block_mat, "rsvd", obj.LR_epsilon); % use the norm to get xi and set precs
                                obj.sq_matrixNorm = obj.sq_matrixNorm + boxK.sq_norm_block_mat{ki};
                            end
                        end
                    end
                end
            end
        end


        function assemble_Dense_Matrices(obj)
            % Dense matrices at leaf level
            j = obj.nLevels + 1;
            nBoxes = obj.nBoxesPerLevel(j);
            for k = 1:nBoxes
                boxK = obj.tree{j}{k};
                if(obj.hodlr_levels == 0 && ~obj.isNBDcompressed)
                    % Neighbors
                    for idx = 1:numel(boxK.neighborNumbers)
                        nn = boxK.neighborNumbers(idx);
                        boxNN = obj.tree{j}{nn};
                        if(obj.is_sym)
                            if(~boxK.isVisited(nn))
                                boxK.dense_nbd_Matrix{nn} = obj.getMatrix(boxK.chargeLocations, boxNN.chargeLocations);
                                if(isempty(boxK.dense_nbd_Matrix{nn}))
                                    continue
                                end
                                norm_dense_nbd_Matrix = norm(boxK.dense_nbd_Matrix{nn}, "fro");
                                obj.sq_matrixNorm = obj.sq_matrixNorm + 2 * (norm_dense_nbd_Matrix * norm_dense_nbd_Matrix);
                                boxNN.isVisited(k) = true;
                            end
                        else
                            boxK.dense_nbd_Matrix{nn} = obj.getMatrix(boxK.chargeLocations, boxNN.chargeLocations);
                            if(isempty(boxK.dense_nbd_Matrix{nn}))
                                continue
                            end
                            norm_dense_nbd_Matrix = norm(boxK.dense_nbd_Matrix{nn}, "fro");
                            obj.sq_matrixNorm = obj.sq_matrixNorm + (norm_dense_nbd_Matrix * norm_dense_nbd_Matrix);
                        end
                    end
                end

                % Self
                boxK.dense_self_Matrix = obj.getMatrix(boxK.chargeLocations, boxK.chargeLocations);
                if(isempty(boxK.dense_self_Matrix))
                    continue;
                end
                norm_dense_self_Matrix = norm(boxK.dense_self_Matrix, "fro");
                obj.sq_matrixNorm = obj.sq_matrixNorm + (norm_dense_self_Matrix * norm_dense_self_Matrix);
            end
        end


        function ap_Hmat_Block_Matrices(obj)
            obj.LR_storage = 0;
            for j = 2:(obj.nLevels-obj.hodlr_levels+1)
                nBoxes = obj.nBoxesPerLevel(j);
                for k = 1:nBoxes
                    boxK = obj.tree{j}{k};

                    % Loop over interacting boxes
                    I_list = boxK.interactionList;
                    for idx = 1:numel(I_list)
                        ki = I_list(idx);
                        if(obj.is_sym)
                            if(~boxK.isVisited(ki))
                                if isempty(boxK.U{ki}) || boxK.sq_norm_block_mat{ki} == 0
                                    continue
                                end
                                xi = sqrt(boxK.sq_norm_block_mat{ki} / obj.sq_matrixNorm); % error
                                update_u = obj.LR_epsilon / (2^((obj.d_dim*(j-1))/2) * xi);

                                find_u = find(obj.unitRoundOff<=update_u);
                                if ~isempty(find_u)
                                    precIndex = obj.sortIdx(find_u(end));
                                else
                                    error('Not available precision');
                                end

                                set_prec(obj.prec_settings{precIndex});
                                [boxK.U{ki}, boxK.V{ki}] = chop_mat(boxK.U{ki}, boxK.S{ki}, boxK.V{ki});  % Chop simulation

                                if size(obj.prec_settings, 2) ~= 0
                                    bits = obj.prec_settings{precIndex}.bits;
                                else
                                    bits = 64;
                                end

                                [m, r] = size(boxK.U{ki});
                                [r, n] = size(boxK.V{ki});

                                obj.LR_storage = obj.LR_storage + (m + n) * r * bits;
                            end
                        else
                            if isempty(boxK.U{ki}) || boxK.sq_norm_block_mat{ki} == 0
                                continue
                            end
                            xi = sqrt(boxK.sq_norm_block_mat{ki} / obj.sq_matrixNorm); % error
                            update_u = obj.LR_epsilon / (2^((obj.d_dim*(j-1))/2) * xi);

                            find_u = find(obj.unitRoundOff<=update_u);

                            if ~isempty(find_u)
                                precIndex = obj.sortIdx(find_u(end));
                            else
                                error('Not available precision');
                            end

                            set_prec(obj.prec_settings{precIndex});
                            [boxK.U{ki}, boxK.V{ki}] = chop_mat(boxK.U{ki}, boxK.S{ki}, boxK.V{ki});  % Chop simulation

                            if size(obj.prec_settings, 2) ~= 0
                                bits = obj.prec_settings{precIndex}.bits;
                            else
                                bits = 64;
                            end

                            [m, r] = size(boxK.U{ki});
                            [r, n] = size(boxK.V{ki});

                            obj.LR_storage = obj.LR_storage + (m + n) * r * bits;

                            obj.Hmat_storage(j) = obj.Hmat_storage(j) + (m + n) * r * bits;
                        end
                    end
                end
                fprintf("LR strorage at level %d is %d\n", j, obj.LR_storage);
            end
        end


        function ap_Junction_Block_Matrices(obj)
            j = obj.nLevels-obj.hodlr_levels+1;
            nBoxes = obj.nBoxesPerLevel(j);
            for k = 1:nBoxes
                boxK = obj.tree{j}{k};
                for idx = 1:numel(boxK.neighborNumbers)
                    nn = boxK.neighborNumbers(idx);

                    if(obj.is_sym)
                        if(~boxK.isVisited(nn))
                            if isempty(boxK.U{nn}) || boxK.sq_norm_block_mat{nn} == 0
                                continue
                            end
                            xi = sqrt(boxK.sq_norm_block_mat{nn} / obj.sq_matrixNorm); % error
                            update_u = obj.LR_epsilon / (2^((obj.d_dim*(j-1))/2) * xi);

                            find_u = find(obj.unitRoundOff<=update_u);
                            if ~isempty(find_u)
                                precIndex = obj.sortIdx(find_u(end));
                            else
                                error('Not available precision');
                            end

                            set_prec(obj.prec_settings{precIndex});
                            [boxK.U{nn}, boxK.V{nn}] = chop_mat(boxK.U{nn}, boxK.S{nn}, boxK.V{nn});  % Chop simulation

                            if size(obj.prec_settings, 2) ~= 0
                                bits = obj.prec_settings{precIndex}.bits;
                            else
                                bits = 64;
                            end

                            [m, r] = size(boxK.U{nn});
                            [r, n] = size(boxK.V{nn});

                            obj.LR_storage = obj.LR_storage + (m + n) * r * bits;
                        end
                    else
                        if isempty(boxK.U{nn}) || boxK.sq_norm_block_mat{nn} == 0
                            continue
                        end
                        xi = sqrt(boxK.sq_norm_block_mat{nn} / obj.sq_matrixNorm); % error
                        update_u = obj.LR_epsilon / (2^((obj.d_dim*(j-1))/2) * xi);

                        find_u = find(obj.unitRoundOff<=update_u);
                        if ~isempty(find_u)
                            precIndex = obj.sortIdx(find_u(end));
                        else
                            error('Not available precision');
                        end

                        set_prec(obj.prec_settings{precIndex});
                        [boxK.U{nn}, boxK.V{nn}] = chop_mat(boxK.U{nn}, boxK.S{nn}, boxK.V{nn});  % Chop simulation

                        if size(obj.prec_settings, 2) ~= 0
                            bits = obj.prec_settings{precIndex}.bits;
                        else
                            bits = 64;
                        end

                        [m, r] = size(boxK.U{nn});
                        [r, n] = size(boxK.V{nn});

                        obj.LR_storage = obj.LR_storage + (m + n) * r * bits;

                        obj.nadm_storage = obj.nadm_storage + (m + n) * r * bits;
                    end
                end
            end
            fprintf("LR strorage at junction level %d is %d\n", j, obj.nadm_storage);
        end


        function ap_HODLR_Block_Matrices(obj)
            for j = (obj.nLevels-obj.hodlr_levels+2):(obj.nLevels+1)
                nBoxes = obj.nBoxesPerLevel(j);
                for k = 1:nBoxes
                    boxK = obj.tree{j}{k};
                    pN = ceil(k / obj.pow2d);
                    nCharges = numel(boxK.charges);
                    % Initialize potential
                    boxK.potential = zeros(nCharges, 1);
                    for c = 1:obj.pow2d
                        ki = obj.tree{j-1}{pN}.childrenNumbers(c);
                        if(ki ~= k)
                            if(obj.is_sym)
                                if(~boxK.isVisited(ki))
                                    if isempty(boxK.U{ki}) || boxK.sq_norm_block_mat{ki} == 0
                                        continue
                                    end
                                    xi = sqrt(boxK.sq_norm_block_mat{ki} / obj.sq_matrixNorm); % error
                                    update_u = obj.LR_epsilon / (2^((obj.d_dim*(j-1))/2) * xi);

                                    find_u = find(obj.unitRoundOff<=update_u);
                                    if ~isempty(find_u)
                                        precIndex = obj.sortIdx(find_u(end));
                                    else
                                        error('Not available precision');
                                    end

                                    set_prec(obj.prec_settings{precIndex});
                                    [boxK.U{ki}, boxK.V{ki}] = chop_mat(boxK.U{ki}, boxK.S{ki}, boxK.V{ki});  % Chop simulation

                                    if size(obj.prec_settings, 2) ~= 0
                                        bits = obj.prec_settings{precIndex}.bits;
                                    else
                                        bits = 64;
                                    end

                                    [m, r] = size(boxK.U{ki});
                                    [r, n] = size(boxK.V{ki});

                                    obj.LR_storage = obj.LR_storage + (m + n) * r * bits;
                                end
                            else
                                if isempty(boxK.U{ki}) || boxK.sq_norm_block_mat{ki} == 0
                                    continue
                                end
                                xi = sqrt(boxK.sq_norm_block_mat{ki} / obj.sq_matrixNorm); % error
                                update_u = obj.LR_epsilon / (2^((obj.d_dim*(j-1))/2) * xi);

                                find_u = find(obj.unitRoundOff<=update_u);
                                if ~isempty(find_u)
                                    precIndex = obj.sortIdx(find_u(end));
                                else
                                    error('Not available precision');
                                end

                                set_prec(obj.prec_settings{precIndex});
                                [boxK.U{ki}, boxK.V{ki}] = chop_mat(boxK.U{ki}, boxK.S{ki}, boxK.V{ki});  % Chop simulation

                                if size(obj.prec_settings, 2) ~= 0
                                    bits = obj.prec_settings{precIndex}.bits;
                                else
                                    bits = 64;
                                end

                                [m, r] = size(boxK.U{ki});
                                [r, n] = size(boxK.V{ki});

                                obj.LR_storage = obj.LR_storage + (m + n) * r * bits;

                                obj.HODLR_strorage(j) = obj.HODLR_strorage(j) + (m + n) * r * bits;
                            end
                        end
                    end
                end
                fprintf("LR strorage at level %d is %d\n", j, obj.HODLR_strorage(j));
            end
        end


        function fast_MVP_Hmat(obj)
            if(obj.N<=obj.N_max)
                obj.Hmat = zeros(obj.N,obj.N); % Initiate
            end
            for j = 2:(obj.nLevels-obj.hodlr_levels+1)
                nBoxes = obj.nBoxesPerLevel(j);
                for k = 1:nBoxes
                    boxK = obj.tree{j}{k};
                    nCharges = numel(boxK.charges);

                    % Initialize potential
                    boxK.potential = zeros(nCharges, 1);

                    % Loop over interacting boxes
                    I_list = boxK.interactionList;
                    for idx = 1:numel(I_list)
                        ki = I_list(idx);
                        boxKI = obj.tree{j}{ki};
                        if(obj.is_sym)
                            if(~boxK.isVisited(ki))
                                if isempty(boxK.U{ki})
                                    continue
                                end

                                set_prec(obj.working_prec);
                                boxK.potential = mchop(boxK.potential + mchop(mchop(boxK.U{ki}) * mchop(mchop(boxK.V{ki}) * mchop(boxKI.charges))));

                                if(obj.N<=obj.N_max)
                                    obj.Hmat(boxK.chargeLocations, boxKI.chargeLocations)  = mchop(boxK.U{ki} * boxK.V{ki});
                                end
                            else
                                if isempty(boxKI.V{k})
                                    continue
                                end

                                set_prec(obj.working_prec);
                                boxK.potential = mchop(boxK.potential + mchop(mchop(boxKI.V{k}') * mchop(mchop(boxKI.U{k}') * mchop(boxKI.charges))));


                                if(obj.N<=obj.N_max)
                                    obj.Hmat(boxK.chargeLocations, boxKI.chargeLocations) = mchop(boxKI.V{k}' * boxKI.U{k}');
                                end
                            end
                        else
                            if isempty(boxK.U{ki})
                                continue
                            end

                            set_prec(obj.working_prec);
                            boxK.potential = mchop(boxK.potential + mchop(mchop(boxK.U{ki}) * mchop(mchop(boxK.V{ki}) * mchop(boxKI.charges))));

                            if(obj.N<=obj.N_max)
                                obj.Hmat(boxK.chargeLocations, boxKI.chargeLocations)  = mchop(boxK.U{ki} * boxK.V{ki});
                            end
                        end
                    end
                end
            end
        end


        function fast_MVP_Junction(obj)
            j = obj.nLevels-obj.hodlr_levels+1;
            nBoxes = obj.nBoxesPerLevel(j);
            for k = 1:nBoxes
                boxK = obj.tree{j}{k};
                for idx = 1:numel(boxK.neighborNumbers)
                    nn = boxK.neighborNumbers(idx);
                    boxNN = obj.tree{j}{nn};

                    if(obj.is_sym)
                        if(~boxK.isVisited(nn))
                            if isempty(boxK.U{nn})
                                continue
                            end

                            set_prec(obj.working_prec);
                            boxK.potential = mchop(boxK.potential + mchop(mchop(boxK.U{nn}) * mchop(mchop(boxK.V{nn}) * mchop(boxNN.charges))));

                            if(obj.N<=obj.N_max)
                                obj.Hmat(boxK.chargeLocations, boxNN.chargeLocations)  = mchop(boxK.U{nn} * boxK.V{nn});
                            end
                        else
                            if isempty(boxNN.V{k})
                                continue
                            end

                            set_prec(obj.working_prec);
                            boxK.potential = mchop(boxK.potential + mchop(mchop(boxNN.V{k}') * mchop(mchop(boxNN.U{k}') * mchop(boxNN.charges))));

                            if(obj.N<=obj.N_max)
                                obj.Hmat(boxK.chargeLocations, boxNN.chargeLocations) = mchop(boxNN.V{k}' * boxNN.U{k}');
                            end
                        end
                    else
                        if isempty(boxK.U{nn})
                            continue
                        end

                        set_prec(obj.working_prec);
                        boxK.potential = mchop(boxK.potential + mchop(mchop(boxK.U{nn}) * mchop(mchop(boxK.V{nn}) * mchop(boxNN.charges))));

                        if(obj.N<=obj.N_max)
                            obj.Hmat(boxK.chargeLocations, boxNN.chargeLocations)  = mchop(boxK.U{nn} * boxK.V{nn});
                        end
                    end
                end
            end
        end


        function fast_MVP_HODLR(obj)
            for j = (obj.nLevels-obj.hodlr_levels+2):(obj.nLevels+1)
                nBoxes = obj.nBoxesPerLevel(j);
                for k = 1:nBoxes
                    boxK = obj.tree{j}{k};
                    pN = ceil(k / obj.pow2d);
                    nCharges = numel(boxK.charges);
                    % Initialize potential
                    boxK.potential = zeros(nCharges, 1);
                    for c = 1:obj.pow2d
                        ki = obj.tree{j-1}{pN}.childrenNumbers(c);
                        if(ki ~= k)
                            boxKI = obj.tree{j}{ki};

                            if(obj.is_sym)
                                if(~boxK.isVisited(ki))
                                    if isempty(boxK.U{ki})
                                        continue
                                    end

                                    set_prec(obj.working_prec);
                                    boxK.potential = mchop(boxK.potential + mchop(mchop(boxK.U{ki}) * mchop(mchop(boxK.V{ki}) * mchop(boxKI.charges))));

                                    if(obj.N<=obj.N_max)
                                        obj.Hmat(boxK.chargeLocations, boxKI.chargeLocations)  = mchop(boxK.U{ki} * boxK.V{ki});
                                    end
                                else
                                    if isempty(boxKI.V{k})
                                        continue
                                    end

                                    set_prec(obj.working_prec);
                                    boxK.potential = mchop(boxK.potential + mchop(mchop(boxKI.V{k}') * mchop(mchop(boxKI.U{k}') * mchop(boxKI.charges))));

                                    if(obj.N<=obj.N_max)
                                        obj.Hmat(boxK.chargeLocations, boxKI.chargeLocations) = mchop(boxKI.V{k}' * boxKI.U{k}');
                                    end
                                end
                            else
                                if isempty(boxK.U{ki})
                                    continue
                                end

                                set_prec(obj.working_prec);
                                boxK.potential = mchop(boxK.potential + mchop(mchop(boxK.U{ki}) * mchop(mchop(boxK.V{ki}) * mchop(boxKI.charges))));

                                if(obj.N<=obj.N_max)
                                    obj.Hmat(boxK.chargeLocations, boxKI.chargeLocations)  = mchop(boxK.U{ki} * boxK.V{ki});
                                end
                            end
                        end
                    end
                end
            end
        end


        function dense_MVP_leaf(obj)
            % Dense MVP at leaf level
            j = obj.nLevels + 1;
            nBoxes = obj.nBoxesPerLevel(j);

            obj.Dense_storage = 0;

            for k = 1:nBoxes
                boxK = obj.tree{j}{k};

                % Neighbors
                if(obj.hodlr_levels == 0 && ~obj.isNBDcompressed)
                    for idx = 1:numel(boxK.neighborNumbers)
                        nn = boxK.neighborNumbers(idx);
                        boxNN = obj.tree{j}{nn};
                        if(obj.is_sym)
                            if(~boxK.isVisited(nn))
                                if(isempty(boxK.dense_nbd_Matrix{nn}))
                                    continue
                                end

                                set_prec(obj.working_prec);
                                boxK.potential = mchop(boxK.potential + mchop(mchop(boxK.dense_nbd_Matrix{nn}) * mchop(boxNN.charges)));

                                obj.Dense_storage = obj.Dense_storage + numel(boxK.dense_nbd_Matrix{nn}) * obj.working_prec.bits;

                                if(obj.N<=obj.N_max)
                                    obj.Hmat(boxK.chargeLocations, boxNN.chargeLocations) = mchop(boxK.dense_nbd_Matrix{nn});
                                end
                            else
                                if(isempty(boxNN.dense_nbd_Matrix{k}))
                                    continue
                                end

                                set_prec(obj.working_prec);
                                boxK.potential = mchop(boxK.potential + mchop(mchop(boxNN.dense_nbd_Matrix{k}') * mchop(boxNN.charges)));

                                if(obj.N<=obj.N_max)
                                    obj.Hmat(boxK.chargeLocations, boxNN.chargeLocations) = mchop(boxNN.dense_nbd_Matrix{k}');
                                end
                            end
                        else
                            if(isempty(boxK.dense_nbd_Matrix{nn}))
                                continue
                            end

                            set_prec(obj.working_prec);
                            boxK.potential = mchop(boxK.potential + mchop(mchop(boxK.dense_nbd_Matrix{nn}) * mchop(boxNN.charges)));

                            obj.Dense_storage = obj.Dense_storage + numel(boxK.dense_nbd_Matrix{nn}) * obj.working_prec.bits;

                            if(obj.N<=obj.N_max)
                                obj.Hmat(boxK.chargeLocations, boxNN.chargeLocations) = mchop(boxK.dense_nbd_Matrix{nn});
                            end
                        end
                    end
                end

                % Self interaction
                if(isempty(boxK.dense_self_Matrix))
                    continue
                end

                set_prec(obj.working_prec);
                boxK.potential = mchop(boxK.potential + mchop(mchop(boxK.dense_self_Matrix) * mchop(boxK.charges)));

                obj.Dense_storage = obj.Dense_storage + numel(boxK.dense_self_Matrix) * obj.working_prec.bits;

                if(obj.N<=obj.N_max)
                    obj.Hmat(boxK.chargeLocations, boxK.chargeLocations) = mchop(boxK.dense_self_Matrix);
                end
            end
        end


        function potential = collectPotential(obj)
            % Collect potentials from all levels
            potential = zeros(obj.N, 1);

            for j = 2:obj.nLevels+1
                startIdx = 1;
                potentialTemp = zeros(obj.N, 1);

                nBoxes = obj.nBoxesPerLevel(j);
                for k = 1:nBoxes
                    boxK = obj.tree{j}{k};
                    len = numel(boxK.potential);
                    potentialTemp(startIdx:startIdx+len-1) = boxK.potential;
                    startIdx = startIdx + len;
                end
                potential = potential + potentialTemp;
            end
        end


        function reorder(obj, potential)
            % Reorder potentials based on charge locations
            potentialTemp = potential;
            startIdx = 1;
            j = obj.nLevels + 1;
            nBoxes = obj.nBoxesPerLevel(j);

            for k = 1:nBoxes
                boxK = obj.tree{j}{k};
                for i = 1:numel(boxK.chargeLocations)
                    idx = boxK.chargeLocations(i);
                    potential(idx) = potentialTemp(startIdx);
                    startIdx = startIdx + 1;
                end
            end

            obj.computed_potential = potential;
        end


        function compute_rel_error_mvp(obj)
            ind = 1:obj.N;
            A = obj.getMatrix(ind,ind);
            true_poten = A * obj.collective_charge;
            errVec = true_poten - obj.computed_potential;
            fprintf("Relative error in MVP:");
            disp(norm(errVec, "fro") / norm(true_poten, "fro"));
        end


        function get_storage(obj)
            fprintf("Low-rank storage\n");
            disp(obj.LR_storage);

            fprintf("Dense matrix storage\n")
            disp(obj.Dense_storage);

            obj.storage = obj.LR_storage + obj.Dense_storage;

            fprintf("Total storage\n")
            disp(obj.storage);
        end


        function compute_Hmat_error(obj)
            H = obj.getMatrix(obj.tree{1}{1}.chargeLocations, obj.tree{1}{1}.chargeLocations);
            nrmH = norm(H, "fro");
            fprintf("Relative error in Matrix construction:\n");
            disp(norm(H-obj.Hmat, "fro")/nrmH);
        end


        function compute_MVP_backward_error(obj)
            H = obj.getMatrix(obj.tree{1}{1}.chargeLocations, obj.tree{1}{1}.chargeLocations);
            nrmH = norm(H, "fro");
            true_mvp = H * obj.collective_charge;
            back_err  = norm(true_mvp - obj.computed_potential, "fro") / (nrmH * norm(obj.collective_charge, "fro"));
            fprintf("Relative backward error in MVP:\n");
            disp(back_err);
        end


        function [sortu, sortIdx] = sort_by_u(~, u_chain)
            callCellFunc = cellfun(@(x)x.u, u_chain);
            [sortu, sortIdx] = sort(callCellFunc);
        end

        %% FOR DEBUGGING PURPOSE
        function visualize_2d_tree(obj)
            figure;
            hold on;
            axis equal;
            axis off;

            x = obj.particle_loc(:,1);
            y = obj.particle_loc(:,2);

            % Plot particles
            scatter(x, y, 10, 'filled', 'b');

            % Plot 2^d-tree boxes recursively with IDs
            plot_box(obj, 1, 1); % start from root: level=1, box=1

            xlabel('x'); ylabel('y');
            hold off;

            function plot_box(obj, level, boxIdx)
                box = obj.tree{level}{boxIdx};
                if isempty(box.chargeLocations)
                    return; % skip empty box
                end

                % box center and half-length
                c = box.center;
                r = obj.boxRadius(level);
                rectangle('Position', [c(1)-r, c(2)-r, 2*r, 2*r], ...
                    'EdgeColor','r','LineWidth',1);

                % Add box ID label
                % text(c(1), c(2), sprintf('%d', box.boxNumber), 'Color', 'k', ...
                %     'HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize',12);

                % recursively plot children
                if level < obj.nLevels+1
                    for cIdx = 1:obj.pow2d
                        childIdx = box.childrenNumbers(cIdx);
                        if childIdx <= numel(obj.tree{level+1})
                            plot_box(obj, level+1, childIdx);
                        end
                    end
                end
            end
        end

        function visualize_3d_tree(obj)
            figure;
            hold on;
            axis equal;
            axis off;
            view(3);

            % Plot particle locations
            scatter3(obj.particle_loc(:,1), obj.particle_loc(:,2), obj.particle_loc(:,3), ...
                10, 'b', 'filled');

            % Loop over all tree levels and boxes
            for j = 1:(obj.nLevels+1)
                for k = 1:obj.nBoxesPerLevel(j)
                    box = obj.tree{j}{k};
                    c = box.center;           % box center
                    r = obj.boxRadius(j);     % box half-length

                    % Create cube vertices
                    [X,Y,Z] = ndgrid([-1 1]*r + c(1), [-1 1]*r + c(2), [-1 1]*r + c(3));
                    verts = [X(:), Y(:), Z(:)];

                    % Define cube edges
                    edges = [1 2; 1 3; 1 5; 2 4; 2 6; 3 4; 3 7; 4 8; 5 6; 5 7; 6 8; 7 8];

                    % Draw edges
                    for e = 1:size(edges,1)
                        plot3([verts(edges(e,1),1) verts(edges(e,2),1)], ...
                            [verts(edges(e,1),2) verts(edges(e,2),2)], ...
                            [verts(edges(e,1),3) verts(edges(e,2),3)], ...
                            'r', 'LineWidth',1.2);
                    end

                    % Label box number at the center
                    % text(c(1), c(2), c(3), sprintf('%d', box.boxNumber), ...
                    %     'Color', 'k', 'HorizontalAlignment','center', ...
                    %     'VerticalAlignment','middle', 'FontSize',12, 'FontWeight','bold');
                end
            end

            view(3);
            grid off;
            axis off;
        end
    end
end
