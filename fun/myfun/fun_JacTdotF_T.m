function [res] = fun_JacTdotF_T(x,FF,theta,tau,Ind,d,idxV,idxF,idxFx,idxFy,...
                            idxFz,Kt,Ae,Ae1x,Ae1y,Ae1z,L,M,N,...
                            circ_L,Aee,K)
%%
num_curr=length(idxF);
Jout=zeros(L,M,N,3);
Jout(idxF)=x(1:num_curr); 
[~,normJ2]=fun_my_postRT2_bis(Jout,Kt,Ae1x,Ae1y,Ae1z,d);
normJ=normJ2.^0.5;
%%
[z_realF]=fun_resistance_NL(normJ,Ind,Kt,L,M,N,Ae,idxF);
z_realx=zeros(L,M,N);
z_realx(idxFx)=z_realF(idxFx);
z_realx_loc=z_realx(idxFx);
z_realy=zeros(L,M,N);
z_realy(idxFy)=z_realF(Kt+idxFy);
z_realy_loc=z_realy(idxFy);
z_realz=zeros(L,M,N);
z_realz(idxFz)=z_realF(2*Kt+idxFz);
z_realz_loc=z_realz(idxFz);
%% 
Nmat=size(Ind,2);
J2_to_dcoeff=zeros(Kt,1);
for ii=1:Nmat
    Ind(ii).ind=reshape(Ind(ii).ind,length(Ind(ii).ind),1);
    if strcmp(Ind(ii).tag,'cond')
        J2_to_dcoeff(Ind(ii).ind,1)=zeros(length(Ind(ii).ind),1);  
    elseif strcmp(Ind(ii).tag,'supercond')
        J2_to_dcoeff(Ind(ii).ind,1)=fun_drho_dJ2(normJ2(Ind(ii).ind),Ind(ii).jc,Ind(ii).ec,Ind(ii).n_exp);
    end
end
[Hjac] = fun_constructH_for_Jacobian_bis(x,d,idxV,idxF,idxFx,idxFy,...
                idxFz,K,Kt,Ae,L,M,N,J2_to_dcoeff,num_curr,normJ2);
%% 
% JAC=[1/tau*L theta*(R+H)       theta*AeeR.';...
%     theta*AeeR          zeros(size(AeeR,1))];
res=multiplyMATVECT_EDDY_NEW5_JAC(FF,theta,tau,...
                                circ_L,z_realx_loc,z_realy_loc,...
                                z_realz_loc,Hjac,...
                                idxFx,idxFy,idxFz,d,Aee,L,M,N);
res=res.';
end

