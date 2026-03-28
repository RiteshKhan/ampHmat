clear
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_dim = 2;

d_dimRootN = [80];

nParticlesInLeafAlong1D = 5;

%LR_epsilon = [1e-2];
LR_epsilon = [1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 1e-1];

kernel_choice = 3;

hodlr_level = 1;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eta = sqrt(d_dim)/2.0;   % eta=1.0 --> HODLRdD (weak), sqrt(d_dim)/2.0 (strong) --> H, 2 HODLR  
particle_location_choice = "uniform"; % uniform, random && 2D: circle, star 3D: sphere, polyhedron. Modifiy in createHmat.m
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


isNBDcompressed = true; % if this flag is flase and hodlr_level = 0 it generates normal H mat.

if(hodlr_level>0)
    isNBDcompressed = true;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
is_sym = false;  % The matrix is symmetric or not
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Working precision
u = precision("d");  % Set working precision

%% %%%%%%%%%%%%%% Mixed precision %%%%%%%%%%%%%%%%%%%%%%
u1 = precision("d");
u2 = precision("s");
u3 = precision("h");
u4 = precision("b");
u5 = precision("q43");
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%% Uniform precision %%%%%%%%%%%%%%%%%%%%%%
% u1 = precision("d");
% u2 = precision("d");
% u3 = precision("d");
% u4 = precision("d");
% u5 = precision("d");
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_chain = prec_chain(u1,u2,u3,u4,u5);

rng(1,'twister');

for i = 1:numel(d_dimRootN)
    for j = 1:numel(LR_epsilon)
        createHmat(d_dim, eta, d_dimRootN(i), hodlr_level, nParticlesInLeafAlong1D, LR_epsilon(j), kernel_choice, particle_location_choice, is_sym, isNBDcompressed, u_chain, u);
    end
end