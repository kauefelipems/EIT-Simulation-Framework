classdef DAC_MODEL
    %Simulates digital-to-analog conversion. Can be used to sample an
    %arbitrary signal vector or directly generate discretized
    %pre-configured signals, such as sine waves.
    %
    %Used to generate PWL signals for PSPICE sources
    
    properties
        fs = 1e6; %sampling frequency
        n_bits = 12; %resolution
        full_scale = 3.3; %full_scale voltage
    end
    
    methods
        
        function obj = DAC_MODEL(fs, n_bits, full_scale)
            obj.fs = fs; 
            obj.n_bits = n_bits;
            obj.full_scale = full_scale;
        end
        
        function out = sample(obj, signal)
            %Samples arbitrary signal
            for i = 1:length(signal)
                sampled_time = signal(i).time(1) : 1/obj.fs : signal(i).time(end);
                
                for j = 1:length(sampled_time)
                    [ d, ix ] = min( abs( signal(i).time - sampled_time(j) ) );
                    out(i).time(j) =  signal(i).time(ix);
                    out(i).amp(j) =  signal(i).amp(ix);
                end
            end                 
        end
        
        function out = discretize(obj, ideal_sample)
            %Discretizes sampled signal amplitude
            for i = 1:length(ideal_sample)                
                offset_sample = obj.full_scale/2 + ideal_sample(i).amp;
                offset_sample = max(min(offset_sample,obj.full_scale),0); %bound values

                LSB = obj.full_scale/(2^obj.n_bits-1); %amplitude discretization
                out(i).amp = LSB*round(offset_sample/LSB);
                out(i).time = ideal_sample(i).time;               
            end         
        end

        function out = sine(obj, f_sig, v_amp, n_periods)
            %Generates DAC sine wave
            
            %D/A time discretization
            DA_time = (0 : 1/obj.fs :(n_periods/f_sig));

            %Signal definition
            ideal_sig = v_amp*sin(2*pi*f_sig*DA_time);
            offset_sig = obj.full_scale/2 + ideal_sig;
            offset_sig = max(min(offset_sig,obj.full_scale),0); %bound values

            %D/A amplitude discretization
            LSB = obj.full_scale/(2^obj.n_bits-1);
            out.amp = transp(LSB*round(offset_sig/LSB));
            out.time = transp(DA_time);
        end
    end
end