%% PREPARATION OF PRECONDITIONER FOR EDDY CURRENT WITH JACOBIAN
% This function has been inspired by 
% lse_sparse_precon_prepare.m which is available at
% https://github.com/acyucel/VoxHenry and in the directory VoxHenry_functions
%
%  [1/tau*L + theta*(R+H)     theta*A'  ] [jc ] 
%  [theta*Ae                      0     ] [phi]  
%%
function [A_inv,LL,UU,PP,QQ,RR,Sch_comp] = preparePREC_EDDY_JAC(theta,tau,d,...
                                        z_realF,idxFx,idxFy,idxFz,...
                                        precL,Hjac,Aee,Kt,prec_type)
%%
LL=[];
UU=[];
PP=[];
QQ=[];
RR=[];
%%
dx = d(1); dy = d(2); dz = d(3);
num_curr = size(Aee,2);
%%
num_currx=length(idxFx);
num_curry=length(idxFy);
num_currz=length(idxFz);
diag_pulse=zeros(num_curr,1);
diag_pulse(1:num_currx,1)                                        =1./(theta*z_realF(idxFx)*dx/(dy*dz)       + precL(1)/tau);
diag_pulse(num_currx+1:num_currx+num_curry,1)                    =1./(theta*z_realF(Kt+idxFy)*dy/(dz*dx)    + precL(2)/tau);
diag_pulse(num_currx+num_curry+1:num_currx+num_curry+num_currz,1)=1./(theta*z_realF(Kt+Kt+idxFz)*dz/(dx*dy) + precL(3)/tau);
inds=zeros(num_curr,3);
inds(1:num_curr,1)=[1:1:num_curr];
inds(1:num_curr,2)=inds(1:num_curr,1);
inds(:,3)=1./(1./(diag_pulse)+diag(theta*Hjac));
A_inv=sparse(inds(:,1),inds(:,2),inds(:,3));
%%
if strcmp(prec_type,'lu')
    Sch_comp=  - (theta*Aee*A_inv*theta*Aee.');
    [LL,UU,PP,QQ,RR] = lu(Sch_comp);
end
end
