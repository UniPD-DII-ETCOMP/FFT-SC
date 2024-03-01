%%
% This function has been inspired by 
% lse_sparse_precon_multiply.m which is available at
% https://github.com/acyucel/VoxHenry
%%
function [x_out] = multiplyPREC(x_in,Ae,A_inv,LL,UU,PP,QQ,RR,prec)
Nn=size(Ae,1); 
Nf=size(Ae,2);
x_out=zeros(Nn+Nf,1);
if strcmp(prec,'lu')
    warning off
    x_out(Nf+1:Nf+Nn) = ...
        QQ * (UU \ (LL \ (PP * (RR \ (x_in(Nf+1:Nf+Nn) - (Ae * (A_inv * x_in(1:Nf))) ) ))));
    warning on
end
x_out(1:Nf) = A_inv*(x_in(1:Nf) - ((Ae.')*x_out(Nf+1:Nf+Nn)));
end