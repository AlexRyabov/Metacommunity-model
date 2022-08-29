function K = MC_GetHalfSaturation(R_eq, r, m)

K = zeros(size(R_eq));
for ri = 1:size(R_eq, 2)
    K(:, ri) = R_eq(:, ri).*(r - m)./m; % half saturation values K
end
end