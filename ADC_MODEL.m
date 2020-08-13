classdef ADC_MODEL
    %Sample data from another sampled dataset with higher sampling 
    %frequency. Used to simulate ADC functionalities on PSPICE signals.
    %
    %Sampling frequency should be at least 100 times smaller than the
    %original sampling. If fs is too close to the original sampling frequency, 
    %the error introduced will be purely cause by the simulation, and
    %therefore results will be incompatible with the experimental measurements.
    %This error can be reduced by interpolation (NOT IMPLEMENTED)
    %
    %PSPICE sampling frequency is not constant, as simulation step changes
    %in time. Therefore, the "maximum simulation step" PSPICE option should 
    %be selected to comply with the desired fs.
    
    properties
        fs = 1e6; %sampling frequency
        n_bits = 16; %resolution
        full_scale = 3.3; %full_scale voltage
    end
    
    methods
        
        function obj = ADC_MODEL(fs, n_bits, full_scale)
            obj.fs = fs; 
            obj.n_bits = n_bits;
            obj.full_scale = full_scale;
        end
        
        function packg = packg(obj, signal, trigger)
            %Creates packages of data to be sampled. Eliminates the data
            %measured outside the sampling window defined by the trigger.
            for i = 1:length(trigger)
                start_pack_index = find(signal.x >= trigger(i).start, 1);
                stop_pack_index  = find(signal.x > trigger(i).stop, 1)-1;
                
                packg(i).time = signal.x(start_pack_index:stop_pack_index);
                packg(i).amp = signal.y(start_pack_index:stop_pack_index);
            end           
        end
        
        function ideal_sample = sample(obj, packg)
            %Sample data packages in time without discretizing the amplitude      
            for i = 1:length(packg)
                sampled_time = packg(i).time(1) : 1/obj.fs : packg(i).time(end);
                
                for j = 1:length(sampled_time)
                    [ d, ix ] = min( abs( packg(1,i).time - sampled_time(j) ) );
                    ideal_sample(i).time(j) =  packg(i).time(ix);
                    ideal_sample(i).amp(j) =  packg(i).amp(ix);
                end
            end                 
        end
        
        function out = discretize(obj, ideal_sample)
            %Amplitude discretization around LSB/2     
            for i = 1:length(ideal_sample)                
                bound_sample = max(min(ideal_sample(i).amp,obj.full_scale),0); %bound values

                LSB = obj.full_scale/(2^obj.n_bits-1); %amplitude discretization
                out(i).amp = LSB*round(bound_sample/LSB);
                out(i).time = ideal_sample(i).time;               
            end
            
        end
        
        function out = avg(obj, dig_signal, n_avg)
            %Averages n_avg peak-to-peak amplitudes
            for i = 1:length(dig_signal)
                out(i) = sum(maxk(dig_signal(i).amp,n_avg)-mink(dig_signal(i).amp,n_avg))./n_avg;
            end
            %making it compatible with inv_solve()
            out = transp(out);
        end

        function out = norm_avg(obj, dig_signal, n_avg)
            %Averages n_avg peak-to-peak amplitudes and normalizes        
            for i = 1:length(dig_signal)
                out(i) = sum(maxk(dig_signal(i).amp,n_avg)-mink(dig_signal(i).amp,n_avg))./n_avg;
            end
            %normalizing and making it compatible with inv_solve()
            out = transp(out/max(out));
        end
    end
end