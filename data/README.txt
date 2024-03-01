-------------------------------------------------------------------------------------------------------
How to create a new user-defined test-case:

1. duplicate the directory test and rename it (e.g. "user_test")
2. copy the user's files "stl1.stl", "stl2.stl",  etc inside "user_test"
3. customize the simulation data between %% BEGIN USER SETTINGS and %% END USER SETTINGS
    Description/Example:
    %% BEGIN USER SETTINGS
    show_mesh = 0; 	          % if 1 make mesh plot
    paraview_export_flag = 1; % if 1 export paraview files in "res_para" directory
    x_ray_flag = 1;           % if 1 make x_ray plot
    model_name='test';         % name of the model
    %
    stl_files(1).name = 'test_cond.stl'; % name of the first stl file to load
    stl_files(1).tag = 'supercond';      % tag for the material (write 'cond' for condutive media, 'supercond' for superconducting media or "terminal" if you want to impose the current)
    stl_files(1).cur=[];                % injected current value, only active if stl_files(1).tag='terminal';
    stl_files(1).rho=1/57e6;            % resistivity of the medium
    
    stl_files(1).name = 'test_cond.stl'; % name of the first stl file to load
    stl_files(1).tag = 'supercond'; % tag for the material (write 'cond' for condutive media, 'supercond' for superconducting media or "terminal" if you want to impose the current)
    stl_files(1).id_e=[]; % id of the related appened element
    stl_files(1).sign=[]; % injected current sign, only active if stl_files(1).tag='terminal';
    stl_files(1).rho=[]; % resistivity of the medium (empty if material is 'supercond')
    stl_files(1).ec=1e-4; % critical field Ec [V/m]
    stl_files(1).n_exp=25; % n exponent E-J power law
    stl_files(1).jc=2.5e6; % critical current Jc [A/m^2]
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
    Nx=25; number of voxels in the x direction
    Ny=10; number of voxels in the y direction
    Nz=10;   number of voxels in the z direction
    % corners
    flag_auto=1; % if 1, user_data below are ignored and the box is created automatically (suggested)
    % user_data   corners of the Box
    meshXmin = 0;      % (m)        
    meshXmax = 1;     % (m)
    meshYmin = -0.015;      % (m)
    meshYmax = 0.015;     % (m)
    meshZmin = -0.015;          % (m)
    meshZmax = 0.015;     % (m)
    %% END USER SETTINGS
3. run "mainVOXELISE.m"

File "data.mat" is now created and you can select "user_dir" in "FFT_SC.m" 
to select the user test case
-------------------------------------------------------------------------------------------------------
