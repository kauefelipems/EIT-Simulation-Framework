%% Setting file path folders

mydir = 'C:\\Users\\Kaue\\Documents\\MATLAB\\EIT Simulation Framework';

netlist_path = [mydir '\\NETLIST files'];
testbench_path = [mydir '\\TESTBENCH files'];
pspice_output_path = [mydir '\\PSPICE files'];

%% Create an SPICE netlist of a FEM model

% These steps, except the eit_spice(), are application dependent.

%Ideal phantom settings
n_elec= 16; 
n_rings= 1;
amp_ideal = 0.5e-3;
med_conductivity = 2e-3;
phantom_conductivity = 0.001e-3;
options = {'no_meas_current','no_rotate_meas'};

%Create EIDORS FEM model (example = 2D circular tank with 16 electrodes)
params= mk_circ_tank(12, [], n_elec );  %EIDORS standard FEM model

params.stimulation= mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}', ...
                            options, amp_ideal); %stimulation structure
params.solve=      'fwd_solve_1st_order';
params.system_mat= 'system_mat_1st_order';

model = eidors_obj('fwd_model', params); %instancing model
show_fem(model ); 

% create homogeneous image + simulated data
mat = med_conductivity * ones( size(model.elems,1) ,1);
homg_img = eidors_obj('image', 'homogeneous image', ...
                     'elem_data', mat, ... 
                     'fwd_model', model); %img_structure
homg_idealdata=fwd_solve(homg_img); %ideal data from the foward solver
                 
% create inhomogeneous image + simulated data
mat([65,81,82,101,102,122])= phantom_conductivity;
inh_img = eidors_obj('image', 'inhomogeneous image', ...
                     'elem_data', mat, ...
                     'fwd_model', model); %img_structure
inh_idealdata=fwd_solve(inh_img); %ideal data from the foward solver

% create PSPICE netlist library (.LIB) at the corresponding folders
% (Basically the eit_spice() function, but generating .lib instead of .s)
eit_pspice(homg_img,[netlist_path '\\homg_net']);       
eit_pspice(inh_img,[netlist_path '\\inhomg_net']); 

%% Set stimulation signal file 

% Create path_list for PWL files
path_list = {};

% Signal parameters
fsignal = 10e3;
v_amp = 1.65;
periods = 5;

% DAC parameters
fs = 1e6;
n_bits = 12;
full_scale = 3.3;

%Generate DAC sampled sine wave 
dac_1 = DAC_MODEL(fs, n_bits, full_scale);
sine_wave = dac_1.sine(fsignal, v_amp, periods);

%Create stimulus file for PWL source
stimulus_file = [testbench_path '\\DA_output.txt'];
path_list = [path_list, stimulus_file];
pwl_write(stimulus_file, sine_wave.time, sine_wave.amp)

%% Set multiplexing patterning files 

%Multiplexer parameters
mux_on = 5;
mux_off = 0;
tsampling = periods/fsignal;
tinj = 1000e-6;
tmeas = 1000e-6;
tinit = 3000e-6;

%Instancing control objects
mux_1 = MUX_CONTROL(mux_on, mux_off, tsampling, tinj, tmeas, tinit);

%Generating control PWL vector and sampling trigger for the ADC
[mux, trigger] = mux_1.pwl_gen(model);

%Construct mux PWL files
for i = 1:ceil(log2(n_elec))
    
    MUX_IP_file = [testbench_path '\\MUX_IP_' int2str(i) '.txt'];
    MUX_IM_file = [testbench_path '\\MUX_IM_' int2str(i) '.txt'];
    MUX_MP_file = [testbench_path '\\MUX_MP_' int2str(i) '.txt'];
    MUX_MM_file = [testbench_path '\\MUX_MM_' int2str(i) '.txt'];

    path_list = [path_list, MUX_IP_file, MUX_IM_file, MUX_MP_file, MUX_MM_file];
    
    pwl_write(MUX_IP_file, mux.time(1:end-1), mux.ip(:,i));
    pwl_write(MUX_IM_file, mux.time(1:end-1), mux.im(:,i));
    pwl_write(MUX_MP_file, mux.time(1:end-1), mux.mp(:,i));
    pwl_write(MUX_MM_file, mux.time(1:end-1), mux.mm(:,i));
end

%% Create file with PWL paths 
% The PWL_paths lists all PWL stimulus paths to facilitate the manual
% source assignment on PSPICE 
PWL_paths = [testbench_path '\\PWL_paths.txt'];
FILE = fopen(PWL_paths, 'wt');
for i=1:length(path_list)
    fprintf(FILE,[path_list{i} '\n']);
end
fclose(FILE);
eidors_msg(['saved PATHS to ' PWL_paths]);

