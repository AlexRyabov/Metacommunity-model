%%
N = 100;
aS = zeros(N, N);

for i = 1:N-1
    aS(i, i + 1) = 1;
    aS(i + 1, i) = 1;
end

    aS(2, end-1) = 0.3;
    aS(end-1, 2) = 0.3;


%aS = aS + 0.05*rand(N, N);
 

aS = aS.^4;

k = 10;
B =(maxk(aS,k));
aSClnd = aS;
aSClnd(aS<repmat(B(end,:), size(aS, 1), 1)) = 0;
aSClnd = max(aSClnd, aSClnd');
aS = aSClnd;
    
L = zeros(size(aS));
for i = 1:size(aS, 2)
    aS(i, i) = 0;
    L(i, :) = -aS(i, :)/sum(aS(i, :), 2);
    %    L(i, :) = -aS(i, :);
end

for i = 1:size(aS, 2)
    L(i, i) = 1;
    %   L(i, i) = sum(aS(i, :), 2);
end


%%
%eigenvalues
[aEV,D] = eig(L);
ev = (diag(D));
[ev, indSortEv] = (sort(ev));
aEV = aEV(:, indSortEv);
%take only first 100 eigenvalues
if length(ev) > 100
    ev = ev(1:100);
    aEV = aEV(:, 1:100);
end

%%
fg01 = f_MakeFigure(1, [1, 1, 1000, 600])
clf
subplot(3, 4, [1,6])
pcolor(aS)
shading flat
cm = gray;
cm = cm(end:-1:1, :);
cm(1, :) = [1, 1, 1];
colormap(cm);
colorbar
title('Similarity')

subplot(3, 4, [3,8])
[~, indSort] = sort(aEV(:, 2))
pcolor(aS(indSort, indSort))
colorbar
shading flat
title('Similarity, sorted')
%%


subplot(3, 4, 9);
plot(ev, '.')
ylabel('Eig val')

subplot(3, 4, 10);
plot(aEV(:, 2), '.')
ylabel('eig vect_2')

subplot(3, 4, 11);
plot(aEV(:, 3), '.')
ylabel('eig vect_3')

subplot(3, 4, 12);
plot(aEV(:, 4), '.')
ylabel('eig vect_4')
