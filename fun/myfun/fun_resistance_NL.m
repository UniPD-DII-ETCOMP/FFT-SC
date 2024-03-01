function [z_realF] = fun_resistance_NL(J_norm,Ind,nVoxel,L,M,N,Ae,idxF)
rhoVoxel=zeros(nVoxel,1);
idxV=[];
Nmat=size(Ind,2);
%%
for ii = 1:Nmat
    Ind(ii).ind=reshape(Ind(ii).ind,length(Ind(ii).ind),1);
    if strcmp(Ind(ii).tag,'cond')
        idxV=[idxV;Ind(ii).ind];  
        rhoVoxel(Ind(ii).ind,1)=Ind(ii).rho;  
    elseif ( strcmp(Ind(ii).tag,'supercond') || strcmp(Ind(ii).tag,'terminal') ) 
        idxV=[idxV;Ind(ii).ind];
        rhoVoxel(Ind(ii).ind,1)=fun_rho_J(J_norm(Ind(ii).ind),Ind(ii).jc,Ind(ii).ec,Ind(ii).n_exp);
    end
end
idxV=unique(idxV);
rho_eV=reshape(rhoVoxel,L,M,N); 
rho_eF=0.5*(abs(Ae(:,:)).'*rho_eV(:)); clear rho_eV
z_realF=rho_eF; clear rho_eF
indFneq=setdiff([1:3*nVoxel].',idxF);
z_realF(indFneq,:)=0;
end