function spice = eit_pspice(img, name)
%function spice = eit_spice(img, [name])
% Converts an EIT FEM model with assigned conductivities (an EIDORS "img") to a
% model reduced, fully connected mesh of resistors in SPICE format.
% If the FEM model has complex valued conductivities, the mesh will be an RLC
% mesh network.
%
% An optional subcircuit 'name' can be provided.
%
% TODO complex value support
% TODO fix electrode ordering for mixed PEM/CEM electrodes
%
% CITATION_REQUEST:
% AUTHOR: A Boyle and A Adler
% TITLE: Integrating Circuit Simulation with EIT FEM Models 
% JOURNAL: 19th International Conference on Biomedical Applications of Electrical Impedance Tomography, Edinburgh, UK
% YEAR: 2018
%

%  (C) 2018 A. Boyle, License: GPL version 2 or version 3

   if ischar(img) & strcmp(img, 'UNIT_TEST') unit_test(); return; end

   if nargin == 1
      name = 'eit';
   end

   Dprime = model_reduce(img);
%  disp(full(1./Dprime))
%  disp(full(-1./(Dprime - tril(Dprime))));
   spice = netlist(Dprime,name);

   if nargout == 0
      filename = [ name '.LIB' ];
      FILE = fopen(filename, 'w');
      fprintf(FILE,'%s\n',spice{:});
      fclose(FILE);
      eidors_msg(['saved SPICE netlist to ' filename]);
      return
   end
end

function Dprime = model_reduce(img)
   Y = calc_system_mat(img); Y= Y.E;
   nn= num_nodes(img);
   % Decompose into blocks; assumes that the nn+1:end nodes are CEM electrodes
   rm = 1:nn; % nodes to fold
   kp = nn+1:size(Y,1); % nodes to keep
   % Now handle PEM electrodes, by transferring nodes between the rm and el sets
   for i = 1:length(img.fwd_model.electrode)
      el = img.fwd_model.electrode(i);
      if length(el.nodes) == 1
         rm(rm == el.nodes) = [];
         kp(end+1) = el.nodes;
      end
   end
   % Note: C = B' ... we don't need to calculate it for symmetric matrices
   A = Y(rm,rm); B= Y(rm,kp); D = Y(kp,kp);
   Dprime  = D - B'*inv(A)*B;
end

function out = netlist(Dprime, name)
   nn = size(Dprime,1);
   ndr = floor(log10(nn*(nn-1)/2))+1; % number of digits for resistors
   nde = floor(log10(nn))+1; % number of digits for electrodes
   str = ['.subckt ' name ];
   for ii = 1:nn;
      str = [ str sprintf([' e%0' num2str(nde) 'd'], ii) ];
   end
   out = { str };
   str = ['re%0' num2str(ndr) 'd e%0' num2str(nde) 'd e%0' num2str(nde) 'd %s'];
   rr = 1;
   for ii = 1:nn;
      for jj = (ii+1):nn;
         val = sprintf('%0.17g',-1/Dprime(ii,jj));
         out(end+1,1) = { strrep(sprintf(str,rr,ii,jj,val),'+','') }; % we strip '+'
         rr = rr +1;
      end
   end
   out(end+1,1) = { '.ends' };
end

function unit_test()
   imdl = mk_common_model('a2s',8);
   stim = mk_stim_patterns(8,1,'{ad}','{ad}',{'meas_current'},1);
   imdl.fwd_model.stimulation = stim(1);
   imdl.fwd_model = rmfield(imdl.fwd_model, 'meas_select');
   img = mk_image(imdl,1);
   v = fwd_solve(img);
   disp('stim');
   disp(full(stim(1).stim_pattern));
   disp(stim(1).stimulation)
   disp('meas');
   disp(full(stim(1).meas_pattern));
   disp('voltages');
   disp(v.meas)
   eit_spice(img)
end
