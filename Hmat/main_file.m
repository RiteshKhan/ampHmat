clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d_dim = 2;

d_dimRootN = 80;

nParticlesInLeafAlong1D = 8;

LR_epsilon = [1e-3 1e-1];

kernel_choice = 2;

eta = sqrt(d_dim)/2.0;

particle_location_choice = "uniform";

is_sym = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tree levels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nLevel = ceil(log(d_dimRootN/nParticlesInLeafAlong1D) / log(2));


% ell = L+1 : standard H-matrix
% ell = L   : hybrid
% ...
% ell = 1   : HODLR

% ell_values = nLevel+1:-1:1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ell_values_array = [nLevel-1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mixed precision
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

useMixedPrecision = true;

u     = precision("d");
u_mvp = precision("d");

if useMixedPrecision
    u1 = precision("d");
    u2 = precision("s");
    u3 = precision("h");
    u4 = precision("b");
    u5 = precision("q52");
else
    u1 = precision("d");
    u2 = precision("d");
    u3 = precision("d");
    u4 = precision("d");
    u5 = precision("d");
end

u_chain = prec_chain(u1,u2,u3,u4,u5);

rng(1,'twister');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileID = fopen('Results/result_amp_hybrid','a+');

fprintf(fileID,'\n\n\n');
fprintf(fileID,'============================================================\n');
fprintf(fileID,'\t\t\t \t\t Start\n');
fprintf(fileID,'============================================================\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Best result MP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

storageHsWP = NaN;

for iN = 1:numel(d_dimRootN)


    % ell = L+1 : standard H-matrix
    % ell = L   : hybrid
    % ...
    % ell = 1   : HODLR

    ell_values = ell_values_array;

    for iEps = 1:numel(LR_epsilon)

        fprintf(fileID,'\n');
        fprintf(fileID,'************* ############################################################ *************\n');
        fprintf(fileID,'d_dimRootN = %d, LR_epsilon = %.2e\n', ...
            d_dimRootN(iN), LR_epsilon(iEps));
        fprintf(fileID,'************* ############################################################ *************\n');

        bestStorageMP = inf;
        bestEllMP = NaN;
        bestHodlrLevelMP = NaN;
        NBDstatus = NaN;

        bestLRStorage = NaN;
        bestDenseStorage = NaN;
        corresWPStorage = NaN;
        bestStorageGainMP = NaN;



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for k = 1:numel(ell_values)

            ell = ell_values(k);

            fprintf('\n');
            fprintf('============================================================\n');
            fprintf('Testing ell = %d / %d\n', ell, nLevel+1);
            fprintf('============================================================\n');

            [totalStorage_MP, LR_storage, Dense_storage, totalStorage_WP, hodlr_level, isNBDcompressed, storage_others] = ...
                createHmat(d_dim, ...
                eta, ...
                d_dimRootN(iN), ...
                ell, ...
                nParticlesInLeafAlong1D, ...
                LR_epsilon(iEps), ...
                kernel_choice, ...
                particle_location_choice, ...
                is_sym, ...
                u_chain, ...
                u, ...
                u_mvp);

            storageGain_MP = totalStorage_WP / totalStorage_MP;

            if(ell == nLevel+1)
                storageHsWP = totalStorage_WP;
            end
            % ---------------------------------------------------------------
            % Save result
            % ---------------------------------------------------------------

            fprintf(fileID, ...
                '\nd_dimRootN = %d, LR_epsilon = %.2e, isMixedPrec = %d, ell = %d, isNBDcompressed = %d, Storage (mixed precision, total) = %.6e GB, Storage (uniform precision, total) = %.6e GB\n', ...
                d_dimRootN(iN), LR_epsilon(iEps), useMixedPrecision, ...
                ell, isNBDcompressed, ...
                totalStorage_MP * 1.25e-10, totalStorage_WP * 1.25e-10);

            fprintf(fileID,"\n*************** == END OF RESULT == ****************\n\n\n");

            % ---------------------------------------------------------------
            % Update best configuration in MP
            % ---------------------------------------------------------------

            if totalStorage_MP < bestStorageMP

                bestStorageMP = totalStorage_MP;

                bestEllMP = ell;
                bestHodlrLevelMP = hodlr_level;
                NBDstatus = isNBDcompressed;

                bestLRStorage = LR_storage;
                bestDenseStorage = Dense_storage;
                corresWPStorage = totalStorage_WP;
                bestStorageGainMP = storageGain_MP;

            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Print optimal configuration in MP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(fileID,'\n');
        fprintf(fileID,'===================================================================================\n');
        fprintf(fileID,'\t\t\t\t OPTIMAL CONFIGURATION (Mixed precision)\n');
        fprintf(fileID,'===================================================================================\n');

        fprintf(fileID,'Optimal ell    ---------------->                                        = %d\n', bestEllMP);
        fprintf(fileID,'HODLR levels   ---------------->                                        = %d\n', bestHodlrLevelMP);
        fprintf(fileID,'NBD compressed                                                          = %d\n', NBDstatus);
        fprintf(fileID,'Low-rank storage                                                        = %.6e GB\n', bestLRStorage * 1.25e-10);
        fprintf(fileID,'Dense storage                                                           = %.6e GB\n', bestDenseStorage * 1.25e-10);
        fprintf(fileID,'Total storage (only payload)                                            = %.6e GB\n', bestStorageMP * 1.25e-10);
        fprintf(fileID,'Corresponding hybrid format Uniform-WP storage                          = %.6e GB\n', corresWPStorage * 1.25e-10);
        fprintf(fileID,'Storage gain (only payload) in hybrid format in MP (over its uniform prec counterpart) = %.6f\n', bestStorageGainMP);
        fprintf(fileID,'Storage gain (only payload) in hybrid format in WP (over uniform prec Hs-matrix)       = %.6f\n', storageHsWP / corresWPStorage);
        fprintf(fileID,'Storage gain (only payload) in hybrid format in MP (over uniform prec Hs-matrix)       = %.6f\n', storageHsWP / bestStorageMP);
        fprintf(fileID,'===================================================================================\n');
    end

end

fclose(fileID);