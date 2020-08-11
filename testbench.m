%% Setting file path folders
pspice_output_path = 'C:\\Users\\Kaue\\Documents\\MATLAB\\EIT Simulation Framework\\PSPICE files\\';
netlist_path = 'C:\\Users\\Kaue\\Documents\\MATLAB\\EIT Simulation Framework\\NETLIST files\\';
testbench_path = 'C:\\Users\\Kaue\\Documents\\MATLAB\\EIT Simulation Framework\\TESTBENCH files\\';

%% Create an SPICE netlist of a FEM model-----------------------------------------------------------------------------------------------

% These steps, except the eit_spice(), are application dependent.

%Ideal phantom settings
n_elec= 16; 
n_rings= 1;
amp_ideal = 0.5e-3;
med_conductivity = 2e-3;
phantom_conductivity = 0.001e-3;
options = {'no_meas_current','no_rotate_meas'};

%Create EIDORS FEM model
params= mk_circ_tank(12, [], n_elec );  

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
% (Basically the eit_spice() function, but generating .LIB instead 
%  of .s files)
eit_pspice(homg_img,[netlist_path 'homg_net']);       
eit_pspice(inh_img,[netlist_path 'inhomg_net']); 

%% Set stimulation signal file -----------------------------------------------------------------------------------------------

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
stimulus_file = [testbench_path 'DA_output.txt'];
path_list = [path_list, stimulus_file];
PWL_write(stimulus_file, sine_wave.time, sine_wave.amp)

%% Set multiplexing patterning files -----------------------------------------------------------------------------------------------

n_measures = length(model.stimulation)*length(model.stimulation(1).meas_pattern(:,1));

%Multiplexers structure
mux.amp = 5;
mux.time = zeros(n_measures+1,1);
mux.ip = zeros(n_measures,ceil(log2(n_elec)));
mux.im = zeros(n_measures,ceil(log2(n_elec)));
mux.mp = zeros(n_measures,ceil(log2(n_elec)));
mux.mm = zeros(n_measures,ceil(log2(n_elec)));

mux.tsampling = periods/fsignal;
mux.tinj = 1000e-6;
mux.tmeas = 1000e-6;
mux.tinit = 3000e-6;
k = 1;
time_step = mux.tinit;

%Construct mux amplitude and time vectors 
for i = 1:length(model.stimulation)
    %Injection and Measurement pattern (1 = I+ and V+ and -1 = I- and V-) 
    inj = model.stimulation(i).stim_pattern/amp_ideal;
    meas = model.stimulation(i).meas_pattern;
    time_step = time_step + mux.tinj;   
    for j = 1:length(model.stimulation(1).meas_pattern(:,1))
        mux.ip(k,:) = mux.amp*de2bi(find(inj==1)-1,ceil(log2(n_elec)));
        mux.im(k,:) = mux.amp*de2bi(find(inj==-1)-1,ceil(log2(n_elec)));
        mux.mp(k,:) = mux.amp*de2bi(find(meas(j,:)==1)-1,ceil(log2(n_elec)));
        mux.mm(k,:) = mux.amp*de2bi(find(meas(j,:)==-1)-1,ceil(log2(n_elec)));
        
        time_step = time_step + mux.tmeas + mux.tsampling;
        mux.time(k+1) = mux.time(k) + time_step;
        trigger(k).start = mux.time(k+1) - mux.tsampling;
        trigger(k).stop = mux.time(k+1);
        time_step = 0;
        k=k+1;        
    end    
end

%Construct mux PWL files
for i = 1:ceil(log2(n_elec))
    
    MUX_IP_file = [testbench_path 'MUX_IP_' int2str(i) '.txt'];
    MUX_IM_file = [testbench_path 'MUX_IM_' int2str(i) '.txt'];
    MUX_MP_file = [testbench_path 'MUX_MP_' int2str(i) '.txt'];
    MUX_MM_file = [testbench_path 'MUX_MM_' int2str(i) '.txt'];

    path_list = [path_list, MUX_IP_file, MUX_IM_file, MUX_MP_file, MUX_MM_file];
    
    PWL_write(MUX_IP_file, mux.time(1:end-1), mux.ip(:,i));
    PWL_write(MUX_IM_file, mux.time(1:end-1), mux.im(:,i));
    PWL_write(MUX_MP_file, mux.time(1:end-1), mux.mp(:,i));
    PWL_write(MUX_MM_file, mux.time(1:end-1), mux.mm(:,i));
end

%% Create file with PWL files paths -----------------------------------------------------------------------------------------------
PWL_paths = [testbench_path 'PWL_paths.txt'];
FILE = fopen(PWL_paths, 'wt');
for i=1:length(path_list)
    fprintf(FILE,[path_list{i} '\n']);
end
fclose(FILE);
eidors_msg(['saved PATHS to ' PWL_paths]);

%% Read the PSPICE output files -----------------------------------------------------------------------------------------------

file_name1 = 'Data-hoimg.csv'
file_name2 = 'Data-inhoimg.csv'
pspice_output.homimg = READ_CURVES([pspice_output_path file_name1]);
pspice_output.inhomimg = READ_CURVES([pspice_output_path file_name2]);

%% Sample and process the measurement vector -----------------------------------------------------------------------------------------------
adc_1 = ADC_MODEL(1e6, 14, 10);

%Digitalize and average homogeneous data
homg_window = adc_1.packg(pspice_output.homimg, trigger);
homg_ideal_samp = adc_1.sample(homg_window);
homg_dig_sample = adc_1.discretize(homg_ideal_samp);
homg_data = adc_1.avg(homg_dig_sample,periods);
homg_data_norm = homg_data/max(homg_data);

%Digitalize and average inhomogeneous data
inh_window = adc_1.packg(pspice_output.inhomimg, trigger);
inh_ideal_samp = adc_1.sample(inh_window);
inh_dig_sample = adc_1.discretize(inh_ideal_samp);
inh_data = adc_1.avg(inh_dig_sample,periods);
inh_data_norm = inh_data/max(inh_data);

%Create structures for data
homg_expdata = homg_idealdata';
homg_expdata.meas = homg_data_norm;
inh_expdata = inh_idealdata';
inh_expdata.meas = inh_data_norm;

%% Reconstruct image -------------------------------------------------------------------------------------------------------------------------

% These steps are application dependent.

% Create model for reconstruction

params= mk_circ_tank(8, [], n_elec ); 

params.stimulation= mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}', ...
                            options, 10);
params.solve=      'fwd_solve_1st_order';
params.system_mat= 'system_mat_1st_order';
params.jacobian=   'jacobian_adjoint';
mdl_2d_2 = eidors_obj('fwd_model', params);
show_fem( mdl_2d_2 );

% Create inverse model
clear inv2d;
inv2d.name= 'EIT inverse';
inv2d.solve=       'np_inv_solve';
inv2d.hyperparameter.value = 3e-3;

inv2d.R_prior= 'prior_TV';
inv2d.reconst_type= 'difference';
inv2d.jacobian_bkgnd.value= 1;
inv2d.fwd_model= mdl_2d_2;
inv2d.fwd_model.misc.perm_sym= '{y}';
inv2d= eidors_obj('inv_model', inv2d);

% Reconstruct and show ideal image
ideal_img= inv_solve( inv2d, inh_idealdata, homg_idealdata);
show_slices(ideal_img);

% Reconstruct and show ideal image
exp_img= inv_solve( inv2d, inh_expdata, homg_expdata);
figure
show_slices(exp_img);