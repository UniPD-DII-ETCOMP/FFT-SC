function [rho] = fun_rho_J(J,Jc,Ec,n_exp)
rho0=1e-2*Ec/Jc;
rhoc=Ec/Jc;
rho = (rhoc*(J/Jc).^(n_exp-1)+rho0);
end