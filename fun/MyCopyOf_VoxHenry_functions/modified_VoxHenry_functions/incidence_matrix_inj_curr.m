%% Incidence Matrix "Lumped current injection"
%**************************************************************************
%%This function computes incidence matrix. 
%%Incidence matrix is relative to faces-to-voxels incidence
%%INPUT: -Kt: Number of voxels in grid (including empty voxels)
%        -LMN: array of number of voxels for each direction (L,M,N)
%        -idxV: index of non-empty voxels
%        -indTerminal_in: index of terminal voxels INPUT
%        -indTerminal_out: index of terminal voxels OUTPUT
%%OUTPUT: -Ae: full incidence matrix
%         -Aee: incidence matrix relative to non-empty voxels
%         -AeeR: reduced incidence matrix to delete rank-deficiency of Aee
%         -idxF: index of faces selected
%         -idxFx: index of faces shared by non-empty volumes along x
%         -idxFy: index of faces shared by non-empty volumes along y
%         -idxFz: index of faces shared by non-empty volumes along z

%%NOTE: Incidence matrix has dimensions: Nvoxels, 3*Nvoxels
%**************************************************************************
function [Ae,Aee,AeeR,idxF,idxFx,idxFy,idxFz,Ae1x,Ae1y,Ae1z] = ...
    incidence_matrix_inj_curr(Kt,LMN,idxV,indTerminal_in,indTerminal_out)
%%Prepare data
L = LMN(1); M = LMN(2); N = LMN(3);

% facciax_[i,j,k] esce dal volume_[i,j,k]
% facciax_[i,j,k] entra nel volume_[i,j,k]+([1,0,0]) 
% facciay_[i,j,k] entra nel volume_[i,j,k]+([0,1,0]) 
% facciaz_[i,j,k] entra nel volume_[i,j,k]+([0,0,1]) 
% seconda della direzione della faccia
Ae=sparse(Kt,3*Kt);
Ae=Ae+sparse(1:Kt,1:Kt,ones(Kt,1),Kt,3*Kt); % facce x shape function ? uscente
Ae=Ae+sparse(1:Kt,Kt+1:2*Kt,ones(Kt,1),Kt,3*Kt); % facce y shape function ? uscente
Ae=Ae+sparse(1:Kt,2*Kt+1:3*Kt,ones(Kt,1),Kt,3*Kt); % facce z shape function ? uscente

idfx_fv=zeros(Kt,2); % volumi attaccati alle facce dove la shape functionx ? entrante
idfy_fv=zeros(Kt,2); % volumi attaccati alle facce dove la shape functiony ? entrante
idfz_fv=zeros(Kt,2); % volumi attaccati alle facce dove la shape functionz ? entrante
for ii = 1:L-1
    for jj = 1:M
        for kk = 1:N
         idfx_fv(ii+(jj-1)*L+(kk-1)*M*L,1)=ii  +(jj  -1)*L+(kk  -1)*M*L;
         idfx_fv(ii+(jj-1)*L+(kk-1)*M*L,2)=ii+1+(jj  -1)*L+(kk  -1)*M*L;
%          idfx_sel(ii+(jj-1)*L+(kk-1)*M*L)=idfx(ii+(jj-1)*L+(kk-1)*M*L)~=0;
        end
    end
end
Ae1x=zeros(2,Kt); 
Ae1x(1,:)=1:Kt;
Ae1x(2,:)=idfx_fv(:,2);

for ii = 1:L
    for jj = 1:M-1
        for kk = 1:N
         idfy_fv(ii+(jj-1)*L+(kk-1)*M*L,1)=ii  +(jj  -1)*L+(kk  -1)*M*L;
         idfy_fv(ii+(jj-1)*L+(kk-1)*M*L,2)=ii  +(jj+1-1)*L+(kk  -1)*M*L;
%          idfy_sel(ii+(jj-1)*L+(kk-1)*M*L)=idfy(ii+(jj-1)*L+(kk-1)*M*L)~=0;
        end
    end
end
Ae1y=zeros(2,Kt);
Ae1y(1,:)=1:Kt;
Ae1y(2,:)=idfy_fv(:,2);

for ii = 1:L
    for jj = 1:M
        for kk = 1:N-1
         idfz_fv(ii+(jj-1)*L+(kk-1)*M*L,1)=ii  +(jj  -1)*L+(kk  -1)*M*L;
         idfz_fv(ii+(jj-1)*L+(kk-1)*M*L,2)=ii  +(jj  -1)*L+(kk+1-1)*M*L;
%          idfz_sel(ii+(jj-1)*L+(kk-1)*M*L)=idfz(ii+(jj-1)*L+(kk-1)*M*L)~=0;
        end
    end
end
Ae1z=zeros(2,Kt);
Ae1z(1,:)=1:Kt;
Ae1z(2,:)=idfz_fv(:,2);

sel=find(idfx_fv(:,2));
Ae=Ae+sparse(idfx_fv(sel,2),idfx_fv(sel,1),-ones(length(sel),1),Kt,3*Kt);
sel=find(idfy_fv(:,2));
Ae=Ae+sparse(idfy_fv(sel,2),Kt+idfy_fv(sel,1),-ones(length(sel),1),Kt,3*Kt);
sel=find(idfz_fv(:,2));
Ae=Ae+sparse(idfz_fv(sel,2),2*Kt+idfz_fv(sel,1),-ones(length(sel),1),Kt,3*Kt);

%%Clear used variables
clear idfx_fv idfx_fv idfz_fv
clear sel


%%Extrapolate incident matrix for non-empty voxels
Aee = Ae(idxV,:); % prendo solo i volumi non vuoti
idxFx = find(sum(abs(Aee(:,1:Kt)),1)==2).'; % facce_x condivise da volumi non vuoti
idxFy = 0*Kt+find(sum(abs(Aee(:,Kt+1:2*Kt)),1)==2).'; % facce_y condivise da volumi non vuoti
idxFz = 0*2*Kt+find(sum(abs(Aee(:,2*Kt+1:3*Kt)),1)==2).'; % facce_z condivise da volumi non vuoti
idxF = find(sum(abs(Aee),1)==2).'; % tutte le facce (x,y,z) considivise dai volumi nn vuoti
% idxF=[idxFx;Kt+idxFy;2*Kt+idxFz];
Aee = Aee(:,idxF); 

%%Incidence matrix is rank deficient, delete at least one row, obtaining "reduced"
%%incidence matrix. This is used in the rest of the code
if isempty(indTerminal_in) && isempty(indTerminal_out)
%     warning('no imposed potential found, last non-empty voxel is
%     grounded') Do not do it when capacitive effects on
    AeeR=Aee;
else
    AeeR_ori=Ae(:,idxF); 
    AeeR=AeeR_ori;
    idxneqV=setdiff([1:Kt].',idxV); % trovo i voxel vuoti
    indTerminal=[indTerminal_in;indTerminal_out]; 
    AeeR_ori(unique([idxneqV;indTerminal]),:)=[]; % butto via i voxel vuoti e butto via i voxel dei terminali
    % correggo la matrice "AeeR2"
    AeeR=[AeeR_ori;sum(AeeR(indTerminal_in,:),1);sum(AeeR(indTerminal_out,:),1)];
end
% tolgo l'ultimo potenziale 
AeeR(end,:)=[]; 

disp(['Size of full incidence matrix: ', num2str(size(Ae))])
disp(['Size of incidence matrix relative to non-empty voxels: ', num2str(size(Aee))])
disp(['Size of reduced incidence matrix: ', num2str(size(AeeR))])

% %% uso la Ae totale per calcolare il rhs dovuto ad eventuali potenziali vincolati
% rhsPOT=Ae(indPotv,idxF).'*valPotv;

end