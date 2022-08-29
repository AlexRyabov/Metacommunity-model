%% Laplace operator 
% This function contains the 2D discreet Laplace operator. Therefore a
% Matrix X is tranformed as described below. Combination of this function
% and the function flux generates flux or no flux boundary conditions. This
% function alone creates periodic boundary conditions. 

function [X_Laplace] = Laplace2D_2x2(X,dx) 
    
[~, ys] = size(X);
if ys == 1
    X_left = [X(:,2:end, :),X(:,1, :)]; % left-shift 
    X_right = [X(:,end, :),X(:,1:(end-1), :)]; % right-shift  
    X_Laplace = ((X_left+X_right-2.*X)./(dx.^2));
else  %2D Laplace operator
    X_up = [X(2:end,:, :);X(1,:, :)]; % up-shift 
    X_left = [X(:,2:end, :),X(:,1, :)]; % left-shift 
    X_down = [X(end,:, :);X(1:(end-1),:, :)]; % down-shift 
    X_right = [X(:,end, :),X(:,1:(end-1), :)]; % right-shift     
    X_Laplace = ((X_left+X_right+X_down+X_up-4.*X)./(dx.^2));
end    
end