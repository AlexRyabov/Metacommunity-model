function S = resourceshuffle(s)
%
%    resourceshuffle
%      This function reshaffles resource supplies over the grid
%      s ia a matrix of resource supplies s(Ly, Lx, Res), where Ly is the
%      vertical size, Lx is the horisontal size and Res is the number of
%      independent resources
%    NOTES
%      
%    HISTORY
%      Alex Ryabov - 2014       : Created.
%      Alex Ryabov - 09/03/2017 : Revised last.
%

 %Reshuffle
 Ly = size(s, 1);
 Lx = size(s, 2);
 Len = Lx*Ly;
 ind = randperm(Len);
 for ri = 1:size(s, 3)
     vS = s(:,:,ri);
     vS = vS(ind);
     S(:,:,ri)=reshape(vS,Ly,Lx);
 end