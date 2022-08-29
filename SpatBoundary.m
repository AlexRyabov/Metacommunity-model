function [X_rr] = SpatBoundary(c,X) 

% This function adds aditional rows and columns around the matrix to get
% different boundary conditions. 
% c = 0, zero flux boundayr conditions
% c = 1 periodical boundary conditions
% Alexey Ryabov,    9.03.209

switch c
    case 0 %% zero flux boundary condition
        X_r = [X(1,:, :); X; X(end,:, :)]; 
        X_rr = [X_r(:,1, :), X_r, X_r(:,end, :)]; 
    case 1 %%periodical boundary conditions
        X_r =  [  X(end,:, :); X;     X(1,:, :)]; 
        X_rr = [X_r(:,end, :), X_r, X_r(:,1, :)]; 

end