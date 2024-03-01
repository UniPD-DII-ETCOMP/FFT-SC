function [fun_drho_dJ2] = fun_drho_dJ2(J2,Jc,Ec,n_exp)
rhoc=Ec/Jc;
const=rhoc/(Jc^(n_exp-1));
fun_drho_dJ2 = const*(n_exp-1)/2*J2.^((n_exp-1)/2-1);
end