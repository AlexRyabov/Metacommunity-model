%Alexey Ryabov 20.03.2019
function [f_N, f_R] = MC_LocalChange(t, R, D, S, N, r, K, m, c, Params)
p = zeros(size(R));  %Limitation of growth by different resources p(x, y, resource#)
gr = zeros(size(N)); %Growth rates in each point gr(x, y, species#)
%Get the number of species and resources
[sp, res] = size(c);  %c(species#, resource#);

%Calculate species growth rates
for si = 1:sp  %loop over species.  c(species#, resource#)
    for ri = 1:res  %loop over resources
        p(:, :, ri) = R(:, :, ri)./(K(si,ri)+ R(:, :, ri)); % Monod functions p
    end
    gr(:, :, si) = r(si)*min(p(:, :, :), [], 3); 
end

f_N = N .* (gr - m);   % change of species biomass

f_R =  (D .* (Params.SSupplActivity(t) * S - R));    % Nutrient change due to nutrient supply

%Calculate consumption
for ri = 1:res
    for si = 1:sp  %loop over species
        f_R(:, :, ri) = f_R(:, :, ri) - (c(si, ri) * (N(:, :, si) .* gr(:, :, si)));
    end;
end;
