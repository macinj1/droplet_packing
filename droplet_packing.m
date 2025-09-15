%% Droplet Packing
% The following code allows to extrac two main information of one figure:
% location (x,y) and the radius of the objects. It is assumeded that the
% figures contains a distribution of circular objects. 
%
% The variable "setting" will save all the parameter needed to reproduce
% the results for a given image. 
%
% The following uses the function rotateAround by Jan Motl (jan@motl.us)
% Jan Motl (2025). Rotate an image around a point (https://www.mathworks.com/matlabcentral/fileexchange/40469-rotate-an-image-around-a-point), MATLAB Central File Exchange. Retrieved September 15, 2025.

clc; clear; close all 

setting.path = '\RawData\20250521_RawImages\' ; 
setting.case = 'P4C3' ; 
setting.file = strcat( setting.case ,  '.tif' ) ;

image = imread( strcat( setting.path , setting.file ) ) ;

%% Rotate image: 
% Select 2 points to rotate the image and place those two points on some
% arbitrary horizontal axis. 

imshow(image)
title('Select 2 points over the SAME edge of the chamber:')
 
[X,Y] = ginput(2);
angle = ( atan( (Y(2) - Y(1)) / (X(2) - X(1)) ) )*180/pi ; 
image = rotateAround( image , X(1), X(2), angle ) ; 
imshow(image)
title('Image rotated')
pause(2)
setting.rotationCenter = [X(1) X(2)];
setting.rotationAngle = angle ;

[I,CroppedSection] = imcrop(image) ;
setting.CroppedSection = CroppedSection ; 
I = rgb2hsv(I) ;
I = I(:,:,3) ;

imshow(I)
title('Cropped Image')
pause(2)
close all 

clear X Y CroppedSection angle 

%% Filtering image
% To filter the type of image used, a maxima intensity is used. To
% determine the best value, a histogram is shown in which the maximum is
% the seed value to later explore. In some cases, it is a good idea to
% sharpen the values. To use it, uncomment the right line. 

J = imcomplement( imlocalbrighten( I ) ) ; 
% setting.sharpenImage = [5 2] ; % [Radius Amount]
% J = imsharpen( imcomplement( imlocalbrighten( I ) ) , Radius=setting.sharpenImage(1) , Amount=setting.sharpenImage(2) ) ; 

histogram(J) 
title('Image intensity values')

%% Using the right thresholding
% The seed value obtained in the previous section is placed in hmin
% parameter. Next, a overlapping image between the original and the one to
% be used in the following sections is shown. The user should explore the
% image and reiterate in case a improvement is needed. 

setting.hmin = 0.202 ; 
BW = logical( watershed( imhmin( J , setting.hmin ) ) ) ; 

imshowpair( I , label2rgb(BW,'abyss','w') , "blend" )  

%% Find objects in the image
% Once the image is processed, four properties are determined: centroid,
% area, circularity and boundary. The last one is used to compute an equivalent 
% circle by fitting the cloud points on the boundary. In addition, droplets
% on the edge are removed to reduce the error in the following
% calculations.

clear area centros cir area_fit radius_fit

BW2 = imclearborder( BW ) ;
stats = regionprops(BW2,'ConvexHull','Area','Centroid','Circularity') ; 

area = cat(1,stats.Area) ; 
cir = cat(1,stats.Circularity) ; 
centros = cat(1,stats.Centroid) ; 

radius_fit = zeros(length(stats),1) ; 
area_fit = zeros(length(stats),1) ; 

for k = 1:length(stats)

    B = stats(k).ConvexHull ; 
    [~,~, radius] = circlefit(B(:,1),B(:,2)) ; 
    radius_fit(k,1) = radius ; 
    area_fit(k,1) = pi*radius^2 ;  
    
end

clear k B radius BW2

%% Filtering objects
% All the objects founded in the image are filtered by using their size and 
% circulary. Modify the following parameter to obtain the best results. 

setting.filter_mean = .9 ;
setting.filter_std  = 0.05 ;
setting.filter_cir  = 0.47 ;
setting.filter_area  = [700 1000000] ;

idx = area_fit./area < setting.filter_mean*mean(area_fit(:,end)./area)+setting.filter_std*std(area_fit(:,end)./area) ... 
    & cir > setting.filter_cir & area > setting.filter_area(1) & area < setting.filter_area(2) ; 

% imshow(BW,[]) ; 
imshow(I)
hold on

% Plotting the circular objects on the image:
plot(centros(idx,1),centros(idx,2),'bo','MarkerFaceColor','b','MarkerSize',3)

%% Droplet definiton
% Create a new variable with the centroids and the radii of the objects
% that will be considered for the analysis. idx should not be deleted! 
% The data obtained is objected centroid, radius and area. 

droplet = cat(1,stats(idx).Centroid) ; 
droplet_area = cat(1,stats(idx).Area) ; 
droplet_radii = sqrt(droplet_area/pi) ; 

clear area area_fit cir radius_fit centros

%% Compute Coefficient of variation (CV)

CV = std(droplet_radii)/mean(droplet_radii) ; 

%% Save data 

clc 
setting.path_save = 'data_processed\' ; 

clearvars -except droplet droplet_radii CV

save( strcat( setting.path_save , setting.case ,  '.mat' ) )
disp( strcat("Case " , setting.case , " was saved." ) )

%% External functions

function output=rotateAround(image, pointY, pointX, angle, varargin)
% ROTATEAROUND rotates an image.
%   ROTATED=ROTATEAROUND(IMAGE, POINTY, POINTX, ANGLE) rotates IMAGE around
%   the point [POINTY, POINTX] by ANGLE degrees. To rotate the image
%   clockwise, specify a negative value for ANGLE.
%
%   ROTATED=ROTATEAROUND(IMAGE, POINTY, POINTX, ANGLE, METHOD) rotates the
%   image with specified method:
%       'nearest'       Nearest-neighbor interpolation
%       'bilinear'      Bilinear interpolation
%       'bicubic'       Bicubic interpolation
%    The default is fast 'nearest'. Switch to 'bicubic' for nicer results.
%
%   Example
%   -------
%       imshow(rotateAround(imread('eight.tif'), 1, 1, 10));
%
%   See also IMROTATE, PADARRAY.
%   Contributed by Jan Motl (jan@motl.us)
%   $Revision: 1.2 $  $Date: 2014/05/01 12:08:01 $
% Parameter checking.
numvarargs = length(varargin);
if numvarargs > 1
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 1 optional input');
end
optargs = {'nearest'};    % Set defaults for optional inputs
optargs(1:numvarargs) = varargin;
[method] = optargs{:};    % Place optional args in memorable variable names
% Initialization.
[imageHeight imageWidth ~] = size(image);
centerX = floor(imageWidth/2+1);
centerY = floor(imageHeight/2+1);
dy = centerY-pointY;
dx = centerX-pointX;
% How much would the "rotate around" point shift if the 
% image was rotated about the image center. 
[theta, rho] = cart2pol(-dx,dy);
[newX, newY] = pol2cart(theta+angle*(pi/180), rho);
shiftX = round(pointX-(centerX+newX));
shiftY = round(pointY-(centerY-newY));
% Pad the image to preserve the whole image during the rotation.
padX = abs(shiftX);
padY = abs(shiftY);
padded = padarray(image, [padY padX]);
% Rotate the image around the center.
rot = imrotate(padded, angle, method, 'crop');
% Crop the image.
output = rot(padY+1-shiftY:end-padY-shiftY, padX+1-shiftX:end-padX-shiftX, :);
end