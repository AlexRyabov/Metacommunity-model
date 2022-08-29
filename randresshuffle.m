function S = randresshuffle (s, r)

 %Reshuffle resources so that supply points with new resource ratios but the same overall
 %amount of resource are produced
 Ly = size(s, 1);
 Lx = size(s, 2);
 Len = Lx*Ly;
 for ri = 1:r
     ind = randperm(Len);
     vS = s(:,:,ri);
     vS = vS(ind);
     S(:,:,ri)=reshape(vS,Ly,Lx);
 end