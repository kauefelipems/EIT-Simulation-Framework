% READ_CURVES CLASS (READ CSV FROM PSPICE)

classdef READ_CURVES
    properties
        y;x
    end
    
    methods
        function obj = READ_CURVES(file)
            curve = csvread(file,2);
            obj.x = curve(:,1);
            obj.y = curve(:,2);
        end
    end
end
