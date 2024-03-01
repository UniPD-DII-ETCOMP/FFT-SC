clear 
clear global
close all
clc
restoredefaultpath
warning on
global meshXmin meshXmax meshYmin meshYmax meshZmin meshZmax
%% BEGIN USER SETTINGS
show_mesh=1;
paraview_export_flag=1;
x_ray_flag=0;
model_name='test';
%
stl_files(1).name='test_cond.stl'; 
stl_files(1).tag='supercond';
stl_files(1).id_e=[]; 
stl_files(1).sign=[];
stl_files(1).rho=[];
stl_files(1).ec=1e-4; 
stl_files(1).n_exp=25; 
stl_files(1).jc=2.5e6;
%
stl_files(2).name='test_terminal1.stl';
stl_files(2).tag='terminal';
stl_files(2).id_e=1;
stl_files(2).sign=1;
stl_files(2).rho=[];
stl_files(2).ec=1e-4;
stl_files(2).n_exp=25;
stl_files(2).jc=2.5e6;
%
stl_files(3).name='test_terminal2.stl';
stl_files(3).tag='terminal';
stl_files(3).id_e=1; 
stl_files(3).sign=-1;
stl_files(3).rho=[];
stl_files(3).ec=1e-4;
stl_files(3).n_exp=25;
stl_files(3).jc=2.5e6;
%
% Box 
% number of voxels in the x y z directions
Nx=25;
Ny=10;
Nz=10;
% corners
flag_auto=1; % if 1, user_data below are ignored
% user_data
meshXmin = 0; % (m)
meshXmax = 1;  % (m)
meshYmin = -0.015;% (m)
meshYmax = 0.015;   % (m)
meshZmin = -0.015;% (m)
meshZmax = 0.015;  % (m)
%% END USER SETTINGS
how_many_stl=size(stl_files,2);
%%
dad=pwd;
cd ..; cd ..; cd('fun'); addpath(genpath(pwd)); cd(dad)
%% Plot the original STL mesh
ccolor=distinguishable_colors(how_many_stl);
figure
hold on
xmin=[];
xmax=[];
ymin=[];
ymax=[];
zmin=[];
zmax=[];
for ii = 1:how_many_stl
[stlcoords] = READ_stl(stl_files(ii).name);
xco = squeeze( stlcoords(:,1,:) )';
yco = squeeze( stlcoords(:,2,:) )';
zco = squeeze( stlcoords(:,3,:) )';
[hpat] = patch(xco,yco,zco,ccolor(ii,:));
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
view(3)
title('stl')
drawnow
% 
xmin=min([xmin;xco(:)]);
xmax=max([xmax;xco(:)]);
ymin=min([ymin;yco(:)]);
ymax=max([ymax;yco(:)]);
zmin=min([zmin;zco(:)]);
zmax=max([zmax;zco(:)]);
%
end
%%
if flag_auto
    meshXmin=xmin;
    meshXmax=xmax;
    meshYmin=ymin;
    meshYmax=ymax;
    meshZmin=zmin;
    meshZmax=zmax;
else
    meshXmin = meshXmin;
    meshXmax = meshXmax; 
    meshYmin = meshYmin; 
    meshYmax = meshYmax;
    meshZmin = meshZmin;
    meshZmax = meshZmax;    
end
%% 
disp('======================================')
disp('Voxelize...')
for ii = 1:how_many_stl
    [o(ii).OUTPUTgrid,...
    o(ii).gridCOx,...
    o(ii).gridCOy,...
    o(ii).gridCOz] = VOXELISE_mod(Nx,Ny,Nz,stl_files(ii).name,'xyz');
o(ii).gridCOx=o(ii).gridCOx;
o(ii).gridCOy=o(ii).gridCOy;
o(ii).gridCOz=o(ii).gridCOz;
end
xyz= grid3dRT2(o(1).gridCOx,o(1).gridCOy,o(1).gridCOz);
xd=xyz(:,:,:,1);
yd=xyz(:,:,:,2);
zd=xyz(:,:,:,3);
%
idx=zeros(Nx,Ny,Nz);
for ii = 1:how_many_stl
idx=idx+(o(ii).OUTPUTgrid);
end
idx=find(idx);
% plot
xidx=xd(idx);
yidx=yd(idx);
zidx=zd(idx);
disp(' ')
%%
disp('======================================')
if x_ray_flag
    disp('X_ray...')
    for ii = 1:how_many_stl
        figure
        subplot(1,3,1);
        title(['xray object',num2str(ii),' ZY'])
        hold on
        imagesc(squeeze(sum(o(ii).OUTPUTgrid,1)));
        colormap(gray(256));
        xlabel('Z-direction');
        ylabel('Y-direction');
        axis equal 
        subplot(1,3,2);
        title(['xray object',num2str(ii),' ZX'])    
        hold on
        imagesc(squeeze(sum(o(ii).OUTPUTgrid,2)));
        colormap(gray(256));
        xlabel('Z-direction');
        ylabel('X-direction');
        axis equal 
        subplot(1,3,3);
        title(['xray object',num2str(ii),' YX'])            
        hold on        
        imagesc(squeeze(sum(o(ii).OUTPUTgrid,3)));
        colormap(gray(256));
        xlabel('Y-direction');
        ylabel('X-direction');
        axis equal 
    end
end
drawnow
%%
if Nx>1
dx=abs(o(1).gridCOx(2)-o(1).gridCOx(1));
else
dx=(max(xco(:))-min(xco(:)));    
end
if Ny>1
dy=abs(o(1).gridCOy(2)-o(1).gridCOy(1));
else
dy=(max(yco(:))-min(yco(:)));    
end
if Nz>1
dz=abs(o(1).gridCOz(2)-o(1).gridCOz(1));
else
dz=(max(zco(:))-min(zco(:)));
end
xidx=xd([idx]);
yidx=yd([idx]);
zidx=zd([idx]);
disp(' ')
%%  paraview
disp('======================================')
if paraview_export_flag
    disp('Paraview...')
P0=[...
    [xidx-dx/2,yidx-dy/2,zidx+dz/2];...
    [xidx-dx/2,yidx-dy/2,zidx-dz/2];...
    [xidx+dx/2,yidx-dy/2,zidx-dz/2];...
    [xidx+dx/2,yidx-dy/2,zidx+dz/2];...
    [xidx-dx/2,yidx+dy/2,zidx+dz/2];...
    [xidx-dx/2,yidx+dy/2,zidx-dz/2];...
    [xidx+dx/2,yidx+dy/2,zidx-dz/2];...
    [xidx+dx/2,yidx+dy/2,zidx+dz/2]];
nVox=length([idx]);
VP=[1:nVox;...
    nVox+1:2*nVox;...
    2*nVox+1:3*nVox;...
    3*nVox+1:4*nVox;...
    4*nVox+1:5*nVox;...
    5*nVox+1:6*nVox;...
    6*nVox+1:7*nVox;...
    7*nVox+1:8*nVox];
val=2*ones(length([idx]),1);
val(1:length(idx))=1;
dad=pwd;
[Ricc] = fun_for_ParaView_sca_HEXA(...
    val,val,P0,VP,dad,model_name);
end
%%
L=length(o(1).gridCOx);
M=length(o(1).gridCOy);
N=length(o(1).gridCOz);
nVoxel=L*M*N;
smeshx=dx;smeshy=dy;smeshz=dz;
%%
for ii = 1:how_many_stl
    Ind(ii).ind = find(o(ii).OUTPUTgrid);
    Ind(ii).tag =stl_files(ii).tag;
    Ind(ii).id_e =stl_files(ii).id_e;
    Ind(ii).sign =stl_files(ii).sign;
    Ind(ii).rho =stl_files(ii).rho;
    Ind(ii).ec = stl_files(ii).ec;
    Ind(ii).n_exp = stl_files(ii).n_exp;
    Ind(ii).jc = stl_files(ii).jc;
end
%%
Nmat = how_many_stl;
save data.mat Ind L M N Nmat nVoxel smeshx smeshy smeshz xyz -v7.3
%%
if show_mesh
   my_show_mesh 
end