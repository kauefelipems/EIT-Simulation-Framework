classdef ADC_MODEL
    
    properties
        fsampling = 1e6; 
        n_bits = 16;
        full_scale = 3.3;
    end
    
    methods
        
        function obj = Set_ADC(obj, fsampling, n_bits, full_scale)
            obj.fsampling = fsampling; 
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
                sampled_time = packg(i).time(1) : 1/obj.fsampling : packg(i).time(end);
                
                for j = 1:length(sampled_time)
                    [ d, ix ] = min( abs( packg(i).time - sampled_time(j) ) );
                    ideal_sample(i).time(j) =  packg(i).time(ix);
                    ideal_sample(i).amp(j) =  packg(i).amp(ix);
                end
            end                 
        end
        
        function dig_signal = digitalize(obj, ideal_sample)

            for i = 1:length(ideal_sample)                
                offset_sample = obj.full_scale/2 + ideal_sample(i).amp;
                offset_sample = max(min(offset_sample,obj.full_scale),0); %bound values

                LSB = obj.full_scale/(2^obj.n_bits-1); %amplitude discretization
                dig_signal(i).amp = LSB*round(offset_sample/LSB) - obj.full_scale/2;
                dig_signal(i).time = ideal_sample(i).time;               
            end
            
        end
        
        function avg_amplitude = avg(obj, dig_signal, n_avg)
            
            for i = 1:length(dig_signal)
                avg_amplitude(i) = sum(maxk(dig_signal(i).amp,n_avg)-mink(dig_signal(i).amp,n_avg))./n_avg;
            end
        end
    end
end