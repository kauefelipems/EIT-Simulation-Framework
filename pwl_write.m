function pwl_write = pwl_write(filename,time,val)
%Write stimulus signal on a .txt file for PWL PSPICE sources
    FILE = fopen(filename, 'wt');
    format = '%.3fus %.3fV\n';
    trise = 1e-9;
    for i = 1:length(time)
        fprintf(FILE,format,time(i)/1e-6,val(i));
        if (i+1 < length(time))
            fprintf(FILE,format,(time(i+1)-trise)/1e-6,val(i));
        end
    end
    fclose(FILE);
    eidors_msg(['saved SPICE PWL to ' filename]);
    return   
end
   
