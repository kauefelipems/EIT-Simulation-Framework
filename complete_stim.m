function stim = complete_stim(n_elec, n_rings, amplitude)
	%Create stimulation structure with all possible electrode combinations
       
	%Setup stimulation pattern
	combin_elec = nchoosek(1:n_elec,2);
	meas_pattern = spalloc(length(combin_elec),n_elec,2);
		
	%Construct stimulation structure
	for i = 1:length(combin_elec)   
		stimulation(i).stimulation = 'Amp';
		stimulation(i).stim_pattern = spalloc(n_elec,1,2);
    
		stimulation(i).stim_pattern(combin_elec(i,1),1) = amplitude ;
		stimulation(i).stim_pattern(combin_elec(i,2),1) = -amplitude;
    
		meas_pattern(i,combin_elec(i,1)) = 1;
		meas_pattern(i,combin_elec(i,2)) = -1;
	end

	for i = 1:length(combin_elec)   
		stimulation(i).meas_pattern = meas_pattern;    
	end

	stim = stimulation;

end