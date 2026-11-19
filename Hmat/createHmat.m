function [totalStorage_MP, LR_storage_MP, Dense_storage_MP, totalStorage_WP, hodlr_level, isNBDcompressed, other_storage] = createHmat(d_dim, eta, d_dimRootN, ell, nParticlesInLeafAlong1D, LR_epsilon, ...
    kernel_choice, ...
    particle_location_choice, ...
    is_sym, ...
    u_chain, ...
    u, ...
    u_mvp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = d_dimRootN^d_dim;

nLevels = ceil(log(d_dimRootN/nParticlesInLeafAlong1D) / log(2));

L = 1;

charges = rand(N,1);

N_max = 10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine format from switching level ell
%
% Convention:
%   ell = 1       -> pure HODLR
%   2 <= ell <= L -> hybrid H/HODLR
%   ell = L+1     -> standard H-matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ell < 1 || ell > nLevels + 1
    error('ell must satisfy 1 <= ell <= nLevels+1.');
end

if ell == 1

    % -------------------------------------------------------------
    % Pure HODLR
    % -------------------------------------------------------------
    hodlr_level      = nLevels;
    isNBDcompressed  = true;

elseif ell == nLevels + 1

    % -------------------------------------------------------------
    % Pure standard H-matrix
    % -------------------------------------------------------------
    hodlr_level      = 0;
    isNBDcompressed  = false;

else

    % -------------------------------------------------------------
    % Hybrid H/HODLR
    %
    % Standard/weak-admissibility part starts at level ell.
    % Finest levels remain HODLR.
    % -------------------------------------------------------------
    hodlr_level      = nLevels - ell;
    isNBDcompressed  = true;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct H-matrix tree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = HmatTree(d_dim, eta, N, N_max, nLevels, hodlr_level, nParticlesInLeafAlong1D, L, LR_epsilon, kernel_choice, is_sym, isNBDcompressed, u_chain, u, u_mvp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tree construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H.createTree();

H.assign_Center_Locations();

H.assign_Tree_Interactions();

%% %%%%%%%%%%%%%% PARTICLE LOCATION (SELECT ONLY ONE) %%%%%%%%%%%%%%%%%%%%
H.assign_Particle_Locations(particle_location_choice); % uniform, random
%H.assign_Particle_Locations_Surface(particle_location_choice); %2D: circle, star 3D: sphere, polyhedron
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Charge locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H.assign_Charge_Locations();

H.assign_NonLeaf_Charge_Locations();

H.assign_Charges(charges);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H-matrix compression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H.compress_Hmat_Block_Matrices();

H.ap_Hmat_Block_Matrices();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Junction compression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isNBDcompressed
    H.compress_Junction_Block_Matrices();
    H.ap_Junction_Block_Matrices();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HODLR compression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if hodlr_level > 0
    H.compress_HODLR_Block_Matrices();
    H.ap_HODLR_Block_Matrices();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dense leaf blocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H.assemble_Dense_Matrices();
H.ap_Dense_NBD_Matrices();

% Matrix-vector product
H.fast_MVP_Hmat();
if(isNBDcompressed)
    H.fast_MVP_Junction(); % Joining level
end
if(hodlr_level>0)
    H.fast_MVP_HODLR();
end

H.dense_MVP_leaf();

poten = H.collectPotential();

H.reorder(poten);

H.calculate_secondary_storage();

H.get_storage();

if(~is_sym)
    fprintf("Hmat storage\n");
    disp(H.Hmat_storage);
    fprintf("Joining level storage\n");
    disp(H.nadm_storage);
    fprintf("HODLR storage\n");
    disp(H.HODLR_strorage);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Return only quantities needed by main file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LR_storage_MP = H.LR_storage;

Dense_storage_MP = H.Dense_storage_MP + H.Dense_storage_diag;

other_storage = H.storage_others;

totalStorage_WP = H.storage_WP + other_storage;

totalStorage_MP = H.storage_MP + other_storage;


fileID = fopen('Results/result_amp_hybrid','a+');

fprintf(fileID,"\n~~~~~~~*************** AMP HMAT ****************~~~~~~~\n");
fprintf(fileID,"Number of particles: %d\n  ", N);
fprintf(fileID,"\nParameters:\n");
fprintf(fileID,"ell value-> %d\t  ", ell);
fprintf(fileID,"LR Tol.-> %d\t  ", LR_epsilon);
fprintf(fileID,"Dim-> %d\t  ", d_dim);
fprintf(fileID,"Min particle at leaf-> %d\t  ", nParticlesInLeafAlong1D^d_dim);
fprintf(fileID,"Kernel choice-> %d\t  ", kernel_choice);
fprintf(fileID,"isNBDcompressed-> %d\t  ", isNBDcompressed);
fprintf(fileID,"isSymmetric-> %d\t  ", is_sym);
fprintf(fileID,"eta-> %d\n", eta);
fprintf(fileID,"Particle location-> %s\n", particle_location_choice);
fprintf(fileID,"=============================\n");

if(N<=N_max)
    txt = evalc('H.compute_rel_error_mvp()');
    fprintf(fileID, "%s", txt);
    fprintf(fileID,"=============================\n");
end

fprintf(fileID,"Level of the tree: ");
fprintf(fileID,"%d\n", nLevels);
fprintf(fileID,"=============================\n");

fprintf(fileID,"Low-rank Storage (GB): ");
fprintf(fileID,"%d\n", LR_storage_MP * 1.25e-10);
fprintf(fileID,"=============================\n");

fprintf(fileID,"Dense matrices Storage (GB): ");
fprintf(fileID,"%d\n", Dense_storage_MP * 1.25e-10);
fprintf(fileID,"=============================\n");

fprintf(fileID,"Total Storage (only payload) (GB): ");
fprintf(fileID,"%d\n", H.storage_MP * 1.25e-10);
fprintf(fileID,"=============================\n");

fprintf(fileID,"Total Storage (total) (GB): ");
fprintf(fileID,"%d\n", (H.storage_MP + H.storage_others) * 1.25e-10);
fprintf(fileID,"=============================\n");

fprintf(fileID,"\n============================= Storage Gain ====================================\n");

fprintf(fileID,"Storage gain (only payload) in MP (over its uniform precision version): %f\n", H.storage_WP/H.storage_MP);
fprintf(fileID,"Storage gain (total) in MP (over its uniform precision version): %f\n", (H.storage_WP + H.storage_others)/(H.storage_MP + H.storage_others));
fprintf(fileID,"===============================================================================\n");


if(N<=N_max)
    txt = evalc('H.compute_Hmat_error()');  % Capture printed output
    fprintf(fileID, "%s\n", txt);
    fprintf(fileID,"=============================\n");
    txt = evalc('H.compute_MVP_backward_error()');  % Backward error in MVP
    fprintf(fileID, "%s", txt);
end

fclose(fileID);   % Close the file

if(~is_sym && d_dim == 2 && N<=8000)
    H.visualize_2d_tree();
elseif(~is_sym && d_dim == 3 && N<=8000)
    H.visualize_3d_tree();
end

% If you want to see the hierarchical matrix
if(~is_sym && N<=8000)
    plotHmatrix(H);
end

end