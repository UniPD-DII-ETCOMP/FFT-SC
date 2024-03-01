%% FFT-SuperConductors
close all
clear global
clear
clc
restoredefaultpath
format short
warning off
delete log_Green.txt
warning on
dad=pwd;
cd('fun'); addpath(genpath(pwd)); cd(dad)
cd('fortran'); addpath(pwd); cd(dad)
%%
%% BEGIN USER SETTINGS
%%
%% Directory
name_dir='test'; % Select folder with data
%% Injected current
Iinj=1500;
%% Frequency
freq=50; 
%% Selections
Integration_flag='NumAn'; %'NumAn'; 'NumNum' (Integration: NumericalNumerical or AnalyticalNumerical)
ext_field_flag=0; % External field
% below you can write the external electric field as a function of x,y,z
% and omega. Active only if ext_field_flag=1
Ex_ext=@(x,y,z,omega) 20e-3*0.5*omega*y;
Ey_ext=@(x,y,z,omega) -20e-3*0.5*omega*x;
Ez_ext=@(x,y,z,omega) 0*z; 
plot_losses_flag=1; % Plot losses over time
%% Solver parameters
fl_precon_type='lu'; % Only LU currently available (Algebraic Multigrid will be available soon)
tol=1e-6;
inner_it=40;
outer_it=5;
NR_tol=1e-5; % Tolerance Newton-Raphson
NR_iter=1000; % Iterations max Newton-Raphson
%% Time integration setting 
tmax=20e-3;
tau=1/freq/250;
time=0:tau:tmax;
theta=1;
%%
%% END USER SETTINGS
%%
%% Load data
cd('data'); cd(name_dir); load('data.mat'); 
fileList = dir('*.stl');
figure
hold on
xmin=[];xmax=[];ymin=[];ymax=[];zmin=[];zmax=[];
ccolor=distinguishable_colors(size(fileList,1));
for ii = 1:size(fileList,1)
    [stlcoords] = READ_stl(fileList(ii).name);
    xco = squeeze( stlcoords(:,1,:) )';
    yco = squeeze( stlcoords(:,2,:) )';
    zco = squeeze( stlcoords(:,3,:) )';
    [hpat] = patch(xco,yco,zco,ccolor(ii,:));
    axis equal
    xlabel('x');
    ylabel('y');
    zlabel('z');
    view(3)
    title('stl (original, not scaled)')
    drawnow
end
cd(dad)
modelname=name_dir;
%%
global flag_equal_dxdydz
%% EM constants
mu0=4*pi*1e-7;
co=299792458;
eps0=1/co^2/mu0;
omega=2*pi*freq(1);
%% Extract data information 
idxV=[]; ind_c=[]; val_c=[]; 
indPotv=[]; valPotv=[];
ape.in=[]; ape.out=[];
for ii = 1:Nmat
    Ind(ii).ind=reshape(Ind(ii).ind,length(Ind(ii).ind),1);
    if strcmp(Ind(ii).tag,'cond')
        idxV=[idxV;Ind(ii).ind];  
    elseif strcmp(Ind(ii).tag,'supercond')
        idxV=[idxV;Ind(ii).ind];
    elseif strcmp(Ind(ii).tag,'terminal')
        idxV=[idxV;Ind(ii).ind]; 
        if Ind(ii).sign==1
            ape(Ind(ii).id_e).in=Ind(ii).ind;
        elseif Ind(ii).sign==-1
            ape(Ind(ii).id_e).out=Ind(ii).ind;
        end
    end
end
idxV=unique(idxV);
%% Grid Definition
disp('======================================')
disp('Grid data...')
disp(['* number of voxels in x direction: ', num2str(L)])
disp(['* number of voxels in y direction: ', num2str(M)])
disp(['* number of voxels in z direction: ', num2str(N)])
disp('* resolutions:')
dx=smeshx; dy=smeshy; dz=smeshz;
if dx==dy && dy==dz
    flag_equal_dxdydz=1;
else
    flag_equal_dxdydz=0;
end
disp([' dx = ',num2str(dx),' m']); 
disp([' dy = ',num2str(dy),' m']);
disp([' dz = ',num2str(dz),' m'])
d=[dx dy dz]; 
Gram=prod(d);
Kt=nVoxel;
K=length(idxV);
%%
disp([' Total number of voxels: ', num2str(Kt)])
disp([' Number of non-empty voxels: ', num2str(K)])
disp(' ')
%% Incidence Matix A
disp('======================================')
disp('Computing incidence matrix...')
mytic_A=tic;
[Ae,Aee,AeeR,idxF,idxFx,idxFy,idxFz,Ae1x,Ae1y,Ae1z] = ...
    incidence_matrix_inj_curr(Kt,[L M N],idxV,ape(1).in,ape(1).out);
mytoc_A=toc(mytic_A);
Nnori=size(Ae,1);
Nn=size(AeeR,1);
Nf=size(AeeR,2);
disp(['Time for computing incidence : ' ,num2str(mytoc_A),' s']);
disp(' ')
%% Compute Green Tensor 
disp('======================================')
disp('Computing Green tensor...')
mytic_G=tic;
[Gmn]=computeGREEN(d,L,M,N,Integration_flag);
mytoc_G=toc(mytic_G);
disp(['Time for computing Green tensor : ' ,num2str(mytoc_G),' s']);
disp(' ')
%% Compute Circulant Tensors
disp('======================================')
disp('Computing circulant tensor...')
mytic_cir=tic;
[circ_L0,preL0]=computeCIRCULANT(Gmn,d,'L'); clear Gmn 
mytoc_circ=toc(mytic_cir);
circ_L0=mu0*circ_L0;
preL0=mu0*preL0;
disp(['Time for computing circulant tensors : ' ,num2str(mytoc_circ),' s'])
disp(' ')
%% Forcing Term: Incident field
disp('======================================')
disp('Computing incident E field...')
mytic_rhs = tic;
if ext_field_flag
    Ex=Ex_ext(xyz(:,:,:,1),xyz(:,:,:,2),xyz(:,:,:,3),omega);
    Ey=Ey_ext(xyz(:,:,:,1),xyz(:,:,:,2),xyz(:,:,:,3),omega);
    Ez=Ez_ext(xyz(:,:,:,1),xyz(:,:,:,2),xyz(:,:,:,3),omega);
    Vx=(Gram.*Ex)./(dy*dz);
    Vy=(Gram.*Ey)./(dx*dz);
    Vz=(Gram.*Ez)./(dx*dy);
    clear Ex Ey Ez     
else
    Vx=zeros(L,M,N);
    Vy=zeros(L,M,N);
    Vz=zeros(L,M,N);  
end
Vext=[Vx(idxFx);Vy(idxFy);Vz(idxFz)];
mytoc_rhs=toc(mytic_rhs);
disp(['Time for computing incident field : ' ,num2str(mytoc_rhs),' s']);
disp(' ')
%% Resistance matrix R(x)
[fun_z_realF]=@(x) fun_resistance_NL(x,Ind,nVoxel,L,M,N,Ae);
%% NL functions
fun.solvJ = @(Xk,FF) fun_resolve_Jac(Xk,FF,theta,tau,Ind,d,idxV,idxF,idxFx,idxFy,...
                            idxFz,Kt,Ae,Ae1x,Ae1y,Ae1z,L,M,N,...
                            circ_L0,AeeR, ...
                            K,preL0,fl_precon_type);
fun.JTFF_T = @(Xk,FF) fun_JacTdotF_T(Xk,FF,theta,tau,Ind,d,idxV,idxF,idxFx,idxFy,...
                            idxFz,Kt,Ae,Ae1x,Ae1y,Ae1z,L,M,N,...
                            circ_L0,AeeR,K);
%% Time-dependent RHS
u=@(tt) [Vext*cos(omega*tt);...
        zeros(Nn-1,1);...
        -Iinj*sin(omega*tt)];
%% Time integration 
disp('====================================')
disp('Time stepping...')
mytic_time=tic;
time=time(:);
ntmax=length(time);
g=zeros(Nf+Nn,1);
gprev=g;
xall=zeros(Nf+Nn,ntmax);
solprev=zeros(Nf+Nn,1);
tt=0;
for nt=1:ntmax
    tt=time(nt);
    disp([' ',num2str(tt),'/',num2str(tmax),' s'])
    g=u(tt);
    RHS = fun_rhs_implicit_theta(solprev,g,gprev,...
                 theta,tau,...
                 Ind,d,idxF,idxFx,idxFy,...
                 idxFz,Kt,L,M,N,...
                 circ_L0,Ae,AeeR,Ae1x,Ae1y,Ae1z);
    fun.F = @(Xk) fun_Ax_minus_b_theta(Xk,solprev,g,gprev,...
                             theta,tau,...
                             Ind,d,idxF,idxFx,idxFy,...
                             idxFz,Kt,L,M,N,...
                             circ_L0,Ae,AeeR,Ae1x,Ae1y,Ae1z);
    options = optimset('TOLFUN',max([NR_tol,NR_tol/10*norm(RHS)]),...
        'MAXITER',NR_iter,'Display','off'); 
    [x,~,~,~,~,~,~,~]=newtonraphson_mod2(fun,solprev*0,options);
    xall(:,nt)=x;
    gprev=g;
    solprev=x;
end
mytoc_time=toc(mytic_time);
disp(['Time for integration : ' ,num2str(mytoc_time),' s'])
disp(' ')
%%
%% POST PROCESSING
%%
%% Losses evaluation
vsol=zeros(Nf+Nn+1,ntmax);
vsol(1:Nf+Nn,:)=xall;
Sout=zeros(ntmax,1);
for nt=1:ntmax
    num_curr=length(idxF);
    Jout=zeros(L,M,N,3);
    Jout(idxF)=vsol(1:Nf,nt); 
    [J,normJ2]=fun_my_postRT2_bis(Jout,Kt,Ae1x,Ae1y,Ae1z,d);
    normJ=normJ2.^0.5;
    rhoVoxel=zeros(nVoxel,1);
    for ii = 1:Nmat
        Ind(ii).ind=reshape(Ind(ii).ind,length(Ind(ii).ind),1);
        if strcmp(Ind(ii).tag,'cond')
            rhoVoxel(Ind(ii).ind,1)=Ind(ii).rho;  
        elseif strcmp(Ind(ii).tag,'supercond') || strcmp(Ind(ii).tag,'terminal')
            rhoVoxel(Ind(ii).ind,1)=fun_rho_J(normJ(Ind(ii).ind),Ind(ii).jc,Ind(ii).ec,Ind(ii).n_exp);
        end
    end
    E=rhoVoxel.*J;
    EdotJ=sum(E.*J,2);
    Sout(nt)=sum(Gram*EdotJ);
end
if plot_losses_flag
    figure
    plot(time,Sout,'-xb','LineWidth',1.5)
    grid on
    xlabel('Time (s)')
    ylabel('AC losses, P_{AC} (W)')
    set(gca,'FontSize',16)
end