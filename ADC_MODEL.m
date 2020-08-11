classdef ADC_MODEL
    
    properties
        fs = 1e6; 
        n_bits = 16;
        full_scale = 3.3;
    end
    
    methods
        
        function obj = ADC_MODEL(fs, n_bits, full_scale)
            obj.fs = fs; 
            obj.n_bits = n_bits;
            obj.full_scale = full_scale;
        end
        
        function packg = packg(obj, signal, trigger)
            
            for i = 1:length(trigger)
                start_pack_index = find(signal.x >= trigger(i).start, 1);
                stop_pack_index  = find(signal.x > trigger(i).stop, 1)-1;
                
                packg(i).time = signal.x(start_pack_index:stop_pack_index);
                packg(i).amp = signal.y(start_pack_index:stop_pack_index);
            end           
        end
        
        function ideal_sample = sample(obj, packg)
            for i = 1:length(packg)
                sampled_time = packg(i).time(1) : 1/obj.fs : packg(i).time(end);
                
                for j = 1:length(sampled_time)
                    [ d, ix ] = min( abs( packg(i).time - sampled_time(j) ) );
                    ideal_sample(i).time(j) =  packg(i).time(ix);
                    ideal_sample(i).amp(j) =  packg(i).amp(ix);
                end
            end                 
        end
        
        function out = discretize(obj, ideal_sample)

            for i = 1:length(ideal_sample)                
                bound_sample = max(min(ideal_sample(i).amp,obj.full_scale),0); %bound values

                LSB = obj.full_scale/(2^obj.n_bits-1); %amplitude discretization
                out(i).amp = LSB*round(bound_sample/LSB);
                out(i).time = ideal_sample(i).time;               
            end
            
        end
        
        function avg_amplitude = avg(obj, dig_signal, n_avg)
            
            for i = 1:length(dig_signal)
                avg_amplitude(i) = sum(maxk(dig_signal(i).amp,n_avg)-mink(dig_signal(i).amp,n_avg))./n_avg;
            end
        end
    end
end