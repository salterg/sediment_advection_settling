function [cf,manning]=compute_friction(h,elev)
friction = zeros(size(h));
hv = zeros(size(h));
m = zeros(size(h));  % Stem density of vegetation [m^-2]

Cb = 65 ; % Background Chezy, default value in Delft3D
CD = 1.2 ; % NOTE: if m*D^2 < 0.006, correction should probably be applied
D = 0.01 ; % Stem diameter [m], assumed
hv=hv+ 0.75 ; % Vegetation height [m], assumed
m((elev >= -0.6) & (elev < -0.04)) = 4;  % Subtidal (Potamogeton, SAVs)
m((elev >= -0.04) & (elev < 0.12)) = 20;  % Low intertidal (Nelumbo, Sagittaria)
m((elev >= 0.12) & (elev < 0.3)) = 50 ; % High intertidal (Colocasia, Polygonum, Scirpus)
m((elev >= 0.3) & (elev < 0.6)) = 100 ; % Supratidal (Salix, high marsh species)

% Baptist-to-Manning's conversion
n_veg = (h.^(1 / 6)).*((Cb.^(-2) + CD .* m .* D .* min(hv,h) ./ 19.62).^(-1 / 2)+ 7.83.* log(max(hv,h)./ hv)).^(-1);
n_veg(h< 0.05) = 0.07; % Enforce rough very shallow water for stability

friction(elev < -10) = 0.02 ; % Wax Lake Outlet, nc
friction(elev >= -10) = 0.015 ; % Bay + WLD Channels, nb
friction(elev >= -0.6) = n_veg(elev >= -0.6) ; % Vegetated areas, subtidal-supratidal
friction(elev >= 0.6) = 0.15 ; % High marsh, above supratidal
        
cf=9.81*(friction.^2).*(h.^(-1/3));
manning=friction;

end