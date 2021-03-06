function [ FluxLim_MINMOD ] = FluxLim_MINMOD( phi_u, phi_c, phi_d )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
DEL = phi_d - phi_u;

if (abs(DEL) < 0.0000001) 
          DWF = 0;  % Upwind
else
          phi_c_norm = (phi_c - phi_u) / DEL;

          if (phi_c_norm < 0)
            DWF = 0; % Upwind
          end
          
          if (phi_c_norm > 1) 
            DWF = 0; % Upwind
          end
          
          if (phi_c_norm >= 0) && (phi_c_norm < 0.5) 
            DWF = 0.5 * phi_c_norm ./ (1 - phi_c_norm);  % Second Order Upwind
          end
          
          if (phi_c_norm >= 0.5) && (phi_c_norm <= 1) 
            DWF = 0.5;  % CD
          end
          
end
       
FluxLim_MINMOD = (1 - DWF) * phi_c + DWF * phi_d;
end
