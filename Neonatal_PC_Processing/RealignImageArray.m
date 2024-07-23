%
% This function converts image array in an arbitrary order into the one
% where the first dimension goes from -x to +x direction, the second
% dimension goes -y to +y direction, and the third dimension goes from -z
% to +z direction.
% 
% In the Philip Achieva system, the x, y, and z directions are the following:
%  from -x to +x : from R to L
%  from -y to +y : from A to P
%  form -z to +z : from F to H
%
% Input:
%    ImageArray : three-dimensional matrix containing image intensity
%    d1 : two-digit character showing the direction of 1st dimension
%    d2 : two-digit character showing the direction of 2nd dimension    
%    d3 : two-digit character showing the direction of 3rd dimension    
%    (d1, d2, and d3 has one of the following: '+x','-x','+y','-y','+z','-z')
%    l_inverse : if -1, perform the inverse conversion
%
% Output:
%    ConvertedArray : three-dimensional matrix of the image where the 1st,
%                     2nd and 3rd dimensions goes x, y, and z directions 
%                     (i.e., from R to L, from A to P, and F to H,
%                     respectively)
%    mat_xyz : 1x3 integer matrix, the sizes of x, y, and z dimensions.
%    xyz_order : 1x3 integer matrix, the dimension of the original image
%                matrix corresponding to x, y, and z respectively.
%    xyz_sign : 1x3 integer matrix, the sign of the original image array 
%               direction w.r.t. x, y, and z.
%
% e.g.) 
% My coronal image from Achieva has the following dimension directions:
%          1st dimension : from R to L --> +x
%          2nd dimension : from H to F --> -z
%          3rd dimension : from A to P --> +y
%
% [ConvertedArray mat_xyz xyr_order] = RealignImageArray(image_array,'+x','-z','+y') 
%
% The output ConvertedArray(mat1,mat3,mat2) gives the image matrix where each dimension 
% goes from R to L (+x), from A to P (+y), and from F to H (+z). 
% mat_xyz = [mat1 mat3 mat2], xyz_order = [1 3 2], xyz_sign = [1 1 -1]
%
% written by Jinsoo Uh (jinsoo.uh@utsouthwestern.edu)
% 2007-02-28

function [f varargout] = RealignImageArray(varargin)
    ImageArray = varargin{1};
    d1 = varargin{2};
    if length(d1) == 6
       dall = d1;
       clear d1;
       d1 = dall(1:2);
       d2 = dall(3:4);
       d3 = dall(5:6);
       if nargin > 2
          l_inverse = varargin{3};
       else
          l_inverse = 1;
       end
    else
       d2 = varargin{3};
       d3 = varargin{4};

       if nargin > 4
          l_inverse = varargin{5};
       else
          l_inverse = 1;
       end
    end
    
    tmp1 = [['1' d1]; ['2' d2]; ['3' d3]];
    tmp2 = sortrows(tmp1,3);
    xyz_order = str2num(tmp2(:,1));
    xyz_sign_pm  = tmp2(:,2);
    
    if l_inverse == 1
       mat1 = size(ImageArray,1);
       mat2 = size(ImageArray,2);
       mat3 = size(ImageArray,3);
       mat123 = [mat1 mat2 mat3];

       mat_x = mat123(xyz_order(1));
       mat_y = mat123(xyz_order(2));
       mat_z = mat123(xyz_order(3));
       mat_xyz = [mat_x mat_y mat_z];

    elseif l_inverse == -1
       mat_x = size(ImageArray,1);
       mat_y = size(ImageArray,2);
       mat_z = size(ImageArray,3);
       mat_xyz = [mat_x mat_y mat_z];

       mat123 = zeros(1,3);
       mat123(xyz_order(1)) = mat_x;
       mat123(xyz_order(2)) = mat_y;
       mat123(xyz_order(3)) = mat_z;
       mat1 = mat123(1);
       mat2 = mat123(2);
       mat3 = mat123(3);              
    end;        
    
    xyz_sign = zeros(1,3);

    if strcmp(xyz_sign_pm(1),'+')
       mat_x_array = 1:mat_x;
       xyz_sign(1) = 1;
    elseif strcmp(xyz_sign_pm(1),'-')            
       mat_x_array = mat_x:-1:1;
       xyz_sign(1) = -1;
    end;

    if strcmp(xyz_sign_pm(2),'+')
       mat_y_array = 1:mat_y;
       xyz_sign(2) = 1;
    elseif strcmp(xyz_sign_pm(2),'-')            
       mat_y_array = mat_y:-1:1;
       xyz_sign(2) = -1;
    end;

    if strcmp(xyz_sign_pm(3),'+')
       mat_z_array = 1:mat_z;
       xyz_sign(3) = 1;       
    elseif strcmp(xyz_sign_pm(3),'-')            
       mat_z_array = mat_z:-1:1;
       xyz_sign(3) = -1;
    end;

    if l_inverse == -1
        ConvertedArray = zeros(mat1,mat2,mat3);
    else
        ConvertedArray = zeros(mat_x,mat_y,mat_z);
    end;
    
    for ii=1:mat1
    for jj=1:mat2
    for kk=1:mat3
        ijk = [ii jj kk];
        if l_inverse == -1
            ConvertedArray(ii,jj,kk)= ImageArray(mat_x_array(ijk(xyz_order(1))), mat_y_array(ijk(xyz_order(2))), mat_z_array(ijk(xyz_order(3))));
        else
            ConvertedArray(mat_x_array(ijk(xyz_order(1))), mat_y_array(ijk(xyz_order(2))), mat_z_array(ijk(xyz_order(3)))) = ImageArray(ii,jj,kk);
        end;
    end;
    end;
    end;
           
    f = ConvertedArray;
    varargout{1} = mat_xyz;
    varargout{2} = xyz_order';
    varargout{3} = xyz_sign;
    
    
    
    
    
    

