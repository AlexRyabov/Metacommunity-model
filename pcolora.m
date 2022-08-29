function pcolora(X, Y, M, XScale, YScale)

%pcolora(X, Y, M, 'log', 'norm')
switch XScale
    case 'log'
        X(end + 1) = X(end) * (X(end) / X(end - 1));
    otherwise
        X(end + 1) = X(end) + (X(end) - X(end - 1));
end

switch YScale
    case 'log'
        Y(end + 1) = Y(end) * (Y(end) / Y(end - 1));
    otherwise
        Y(end + 1) = Y(end) + (Y(end) - Y(end - 1));
end
% add dummy column and dummy row
M = [M, NaN(size(M, 1), 1)];
M = [M; NaN(1, size(M, 2))];
% M = [NaN(1, size(M, 2));M];
% M = [NaN(size(M, 1), 1), M];

pcolor(X, Y, M)
set(gca, 'XScale', XScale, 'YScale', YScale);

shading flat
end