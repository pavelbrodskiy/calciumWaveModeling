%[Xq,Yq] = meshgrid([0.1:0.1:0.5],[1:1:7]);
%Vq = interp2(X,log(Y),V,Xq,Yq);

% [x,y] = meshgrid([0.1:0.1:0.5],[1:1:7]);
% tri = delaunay(x,y);
% z = peaks(15);
% trisurf(tri,x,y,z)

y2 = Y;
y2(y2>100) = 100;
scatter(X,y2,20,p);