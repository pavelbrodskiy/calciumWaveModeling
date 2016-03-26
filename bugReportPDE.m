% This is to demonstrate a bug in function checkDCoefSize and possibly the
% other size-checking functions.
%
% Error occurs when testing the size of d in CoefficientAssignment.m
% 
% function checkDCoefSize(dcoef, systemsize)
%     if isscalar(dcoef) && ~isa(dcoef, 'function_handle') && dcoef == 0
%         return
%     elseif isa(dcoef, 'function_handle')
%         return
%     end
%     dveclen = numel(dcoef);
%     lengthok = (dveclen == 0 || dveclen == 1 || dveclen == systemsize || ...
%         dveclen == (systemsize*(systemsize+1)/2) || dveclen == systemsize*systemsize);
%     if ~(isvector(dcoef)  && iscolumn(dcoef) && lengthok)
%         error(message('pde:pdeCoefficientSpecification:invalidDMatrixSize'));
%     end
% end
% 
% Appears to be because isvector(dcoef) && iscolumn(dcoef) evaluates to 0
% in all cases other than if d is a scalar. However, lengthok considers
% systemsize*systemssize (such as an NxN matrix) to be a suitable size. And
% lengthok is redundant if this is meant to check to see if d is a matrix.
%
% isvector() && iscolumn() is repeated in numerous locations, meaning that
% m, d, c, a, and f must all be scalar inputs despite the manual listing
% them as matrix inputs.

clearvars

% Initialize system of 4 PDEs on a square domain of size 1
model = createpde(4);
dl = decsg([3;4;0;1;1;0;0;0;1;1],'SQ1',[83;81;49]);
geometryFromEdges(model, dl);

% Apply no-flux boundary conditions
applyBoundaryCondition(model,'edge',1:model.Geometry.NumEdges,'g',0,'q',0);

% Set initial conditions to 1
setInitialConditions(model, [1; 1; 1; 1]);

% Set trivial system coefficients
m = 0;
d = eye(4);
c = eye(4);
a = eye(4);
f = 0;

specifyCoefficients(model,'m',m,'d',d,'c',c,'a',a,'f',f);

% Solve system
generateMesh(model);
results = solvepde(model, ts);
u = results.NodalSolution;
pdeplot3D(model,'colormapdata',u)