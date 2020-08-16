%% Read the PSPICE output files 
file_name1 = '\\Data-hoimg.csv';
file_name2 = '\\Data-inhoimg.csv';
pspice_output.homimg = READ_CURVES([pspice_output_path file_name1]);
pspice_output.inhomimg = READ_CURVES([pspice_output_path file_name2]);

%% Sample and process the measurement vector 
adc_1 = ADC_MODEL(1e6, 14, 10);

%Digitalize and average homogeneous data
homg_window = adc_1.packg(pspice_output.homimg, trigger);
homg_ideal_samp = adc_1.sample(homg_window);
homg_dig_sample = adc_1.discretize(homg_ideal_samp);
homg_data_norm = adc_1.norm_avg(homg_dig_sample,periods);

%Digitalize and average inhomogeneous data
inh_window = adc_1.packg(pspice_output.inhomimg, trigger);
inh_ideal_samp = adc_1.sample(inh_window);
inh_dig_sample = adc_1.discretize(inh_ideal_samp);
inh_data_norm = adc_1.norm_avg(inh_dig_sample,periods);


%% Reconstruct image 

% Create EIDORS measurement structures for PSPICE data
homg_expdata.meas = homg_data_norm;
homg_expdata.time = NaN;
homg_expdata.name = 'solved by fwd_solve_1st_order';
homg_expdata.type = 'data';

inh_expdata.meas = inh_data_norm;
inh_expdata.time = NaN;
inh_expdata.name = 'solved by fwd_solve_1st_order';
inh_expdata.type = 'data';

%These steps are application dependent.

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

% Reconstruct and show experimental image
exp_img= inv_solve( inv2d, inh_expdata, homg_expdata);
figure
show_slices(exp_img);
