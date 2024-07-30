function map = D2GaussFunctionRot(parameters, mesh)
%D2GAUSSFUNCTIONROT   Create 2D Gaussian map using the given parameters and
% size of mesh.

% INPUTS
% parameters    [amplitude, x-center, width, y-center, height, rotation]
% mesh          [rows x columns x 2], mesh(:,:,1): x-coordinates,
%               mesh(:,:,2): y-coordinates

% OUTPUTS
% map           [rows x columns], map of 2D Gaussian

% rotate input matrix of x- and y-coordinates
xdatarot(:,:,1)= mesh(:,:,1)*cos(parameters(6)) - mesh(:,:,2)*sin(parameters(6));
xdatarot(:,:,2)= mesh(:,:,1)*sin(parameters(6)) + mesh(:,:,2)*cos(parameters(6));
% coordinates pointing to center of Gaussian
x0rot = parameters(2)*cos(parameters(6)) - parameters(4)*sin(parameters(6));
y0rot = parameters(2)*sin(parameters(6)) + parameters(4)*cos(parameters(6));

map = parameters(1) * exp( ...
    -((xdatarot(:,:,1) - x0rot).^2 / (2 * parameters(3)^2) + ...
      (xdatarot(:,:,2) - y0rot).^2 / (2 * parameters(5)^2)));
