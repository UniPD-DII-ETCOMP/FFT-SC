%% System MVM with vector: Eddy Currents Problem theta-tau
%**************************************************************************
%%This function computes MVM of system matrix and a vector
%System matrix is:
%                       [R+iwL Aee']
%                       [Aee    0  ]
%Ref: [1] "VoxHenry: FFT-Accelerated Inductance Extraction for Voxelized
%Geometries"
%Ref: [2] "VoxHenry-master"
%%INPUT: -JIn: (Current+Potential) vector with dimension equal to selected
%%faces & selected volume/nodes
%        -fN: FFT-Circulant Tensor of relative to L matrix
%        - z_real: resistance matrix (x,y,z)
%        - idxFx,y,z: FACE index of non-empty voxels 3 directions
%        - d = (dx,dy,dz): array of voxel dimensions
%        - AeeR: REDUCED incidence matrix for selected faces, volumes
%        -L,M,N total number of voxels for each direction
%%OUTPUT: - JOut: (Current+Potential) output vector
%**************************************************************************
function [sol] = multiplyMATVECT_EDDY_THTAU(JIn0,theta,tau,fN,z_realx_loc,z_realy_loc,...
                                    z_realz_loc,idxFx,idxFy,idxFz,d,AeeR,L,M,N)
global flag_equal_dxdydz

%%Obtain grid steps
dx = d(1); dy = d(2); dz = d(3);
% Kt = L*M*N; %total number of voxels

%%Dimensions of Problem
num_node=size(AeeR,1); %num potential nodes
num_curr=size(AeeR,2); %num current faces

%%Allocate space
% JIn = zeros(L, M, N, 3); %3 because we have 3 basis functions
% Jout = zeros(2*L, 2*M, 2*N);    
% translate from local (idx) to global (L,M,N) coordinates
% JIn([idxFx;idxFy+Kt;idxFz+2*Kt]) = JIn0(1:num_curr);

sol=zeros(num_curr+num_node,1);

%%1.*************************** Perform Z*j ***************************
%%Apply fft and mv-op for each of the components of JIn
%x component of JIn, store contribution on 3 components of Jout
if dx == dy && dy == dz && flag_equal_dxdydz
    indloc = 1:length(idxFx);
    fJ=zeros(L, M, N);
    fJ(idxFx) = JIn0(indloc);
    fJ = fftn(fJ,[2*L, 2*M, 2*N]);
    Jout = (fN(:,:,:,1)) .* fJ; % Gxx*Jx
    clear fJ
    Jout = ifftn(Jout);
    Jout =Jout(1:L,1:M,1:N);
    sol(indloc)=Jout(idxFx);
    clear Jout
    sol(indloc) = theta * (dx/(dy*dz)) * (z_realx_loc .* JIn0(indloc)) + sol(indloc)/tau;
  
    indloc = length(idxFx)+1:length(idxFx)+length(idxFy);
    fJ=zeros(L, M, N);
    fJ(idxFy) = JIn0(indloc);
    fJ = fftn(fJ,[2*L, 2*M, 2*N]);
    Jout = (fN(:,:,:,1)) .* fJ; % Gxx*Jx
    clear fJ
    Jout = ifftn(Jout);
    Jout =Jout(1:L,1:M,1:N);
    sol(indloc)=Jout(idxFy);
    clear Jout
    sol(indloc) = theta * (dy/(dx*dz)) * (z_realy_loc .* JIn0(indloc)) + sol(indloc)/tau;

    indloc = length(idxFx)+length(idxFy)+1:length(idxFx)+length(idxFy)+length(idxFz);
    fJ=zeros(L, M, N);
    fJ(idxFz) = JIn0(indloc);
    fJ = fftn(fJ,[2*L, 2*M, 2*N]);
    Jout = (fN(:,:,:,1)) .* fJ; % Gxx*Jx
    clear fJ
    Jout = ifftn(Jout);
    Jout =Jout(1:L,1:M,1:N);
    sol(indloc)=Jout(idxFz);
    clear Jout
    sol(indloc) = theta * (dz/(dx*dy)) * (z_realz_loc .* JIn0(indloc)) + sol(indloc)/tau;
else
    indloc = 1:length(idxFx);
    fJ=zeros(L, M, N);
    fJ(idxFx) = JIn0(indloc);
    fJ = fftn(fJ,[2*L, 2*M, 2*N]);
    Jout = (fN(:,:,:,1)) .* fJ; % Gxx*Jx
    clear fJ
    Jout = ifftn(Jout);
    Jout =Jout(1:L,1:M,1:N);
    sol(indloc)=Jout(idxFx);
    clear Jout
    sol(indloc) = theta * (dx/(dy*dz)) * (z_realx_loc .* JIn0(indloc)) + sol(indloc)/tau;
  
    indloc = length(idxFx)+1:length(idxFx)+length(idxFy);
    fJ=zeros(L, M, N);
    fJ(idxFy) = JIn0(indloc);
    fJ = fftn(fJ,[2*L, 2*M, 2*N]);
    Jout = (fN(:,:,:,2)) .* fJ; % Gxx*Jx
    clear fJ
    Jout = ifftn(Jout);
    Jout =Jout(1:L,1:M,1:N);
    sol(indloc)=Jout(idxFy);
    clear Jout
    sol(indloc) = theta * (dy/(dx*dz)) * (z_realy_loc .* JIn0(indloc)) + sol(indloc)/tau;

    indloc = length(idxFx)+length(idxFy)+1:length(idxFx)+length(idxFy)+length(idxFz);
    fJ=zeros(L, M, N);
    fJ(idxFz) = JIn0(indloc);
    fJ = fftn(fJ,[2*L, 2*M, 2*N]);
    Jout = (fN(:,:,:,3)) .* fJ; % Gxx*Jx
    clear fJ
    Jout = ifftn(Jout);
    Jout =Jout(1:L,1:M,1:N);
    sol(indloc)=Jout(idxFz);
    clear Jout
    sol(indloc) = theta * (dz/(dx*dy)) * (z_realz_loc .* JIn0(indloc)) + sol(indloc)/tau;    
end

%%2.*************************** Perform A'*phi ***************************
sol(1:num_curr) = sol(1:num_curr) + theta * (AeeR.'*JIn0(num_curr+1:num_curr+num_node)) ;

%%3.*************************** Perform A*j ***************************
sol(num_curr+1:num_curr+num_node) = theta * AeeR*JIn0(1:num_curr);

end

