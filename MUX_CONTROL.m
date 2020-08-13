classdef MUX_CONTROL
    %Generation and management of stimulus signals for the switching 
    %control bus
    
    properties
        v_on = 5; % control voltage value for logic 1
        v_off = 0; % control voltage value for logic 0
        t_samp = 0.5e-3; % sampling time (used to generate the 
                         % ADC trigger signal)
        t_inj = 1e-3; % stabilization time for injection switching 
        t_meas = 1e-3; % stabilization time for measurement switching
        t_init = 1e-3; % initialization time delay
    end
    
    methods
        
        function obj = MUX_CONTROL(v_on, v_off,...
                                        t_samp, t_inj, t_meas, t_init)
            obj.v_on = v_on;
            obj.v_off = v_off; 
            obj.t_samp = t_samp; 
            obj.t_inj = t_inj; 
            obj.t_meas = t_meas;
            obj.t_init = t_init; 
            
        end
        
        function [mux, trigger] = pwl_gen(obj, model)
            % Creates stimulus vectors for the PWL PSPICE source,
            % Also generates a trigger file to be used by the ADC_MODEL
            % class for packaging the signal

            n_injections = length(model.stimulation);
            n_voltages = length(model.stimulation(1).meas_pattern(:,1));
            n_meas = n_injections * n_voltages;
            n_elec = length(model.electrode);
            
            % Initialize output mux structure
            % ip = injection +, im = injection -
            % mp = voltage +, mm = voltage -
            mux.time = zeros(n_meas+1,1);
            mux.ip = zeros(n_meas,ceil(log2(n_elec)));
            mux.im = zeros(n_meas,ceil(log2(n_elec)));
            mux.mp = zeros(n_meas,ceil(log2(n_elec)));
            mux.mm = zeros(n_meas,ceil(log2(n_elec)));
            
            k = 1;
            time_step = obj.t_init;

            %Sweep through injection electrodes sequence
            for i = 1 : n_injections
                
                %Get current and measurement electrodes for each injection
                inj = model.stimulation(i).stim_pattern;
                meas = model.stimulation(i).meas_pattern;
                time_step = time_step + obj.t_inj;   
                
                %Sweep through measurement electrodes sequence
                for j = 1 : n_voltages
                    
                    %Build control voltage and time vectors (each mux control
                    %is a bus with ceil(log2(n_elec)) bit signals, for n_elec
                    %binary representation)
                    mux.ip(k,:) = obj.v_on*de2bi(find(inj > 0)-1,ceil(log2(n_elec)));
                    mux.im(k,:) = obj.v_on*de2bi(find(inj < 0)-1,ceil(log2(n_elec)));
                    mux.mp(k,:) = obj.v_on*de2bi(find(meas(j,:) > 0)-1,ceil(log2(n_elec)));
                    mux.mm(k,:) = obj.v_on*de2bi(find(meas(j,:) < 0)-1,ceil(log2(n_elec)));
        
                    time_step = time_step + obj.t_meas + obj.t_samp;
                    mux.time(k+1) = mux.time(k) + time_step;
                    
                    %Build trigger signal for ADC sampling
                    trigger(k).start = mux.time(k+1) - obj.t_samp;
                    trigger(k).stop = mux.time(k+1);
                    
                    time_step = 0;
                    k=k+1;        
                end    
                
                %Replace zeros on the binary bus with the v_off voltage          
                mux.ip(mux.ip == 0) = obj.v_off;
                mux.im(mux.im == 0) = obj.v_off;
                mux.mp(mux.mp == 0) = obj.v_off;
                mux.mm(mux.mm == 0) = obj.v_off;                
            end
            
        end
    end
end