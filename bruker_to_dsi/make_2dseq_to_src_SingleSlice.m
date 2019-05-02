file = '/Volumes/Alan_RAID/TCD/LallyLabDesktop/Data/20190225_BT/18/pdata/1/2dseq';

num_dif = 11;
dimension = [96 96 1]; % [x y z] the dimension of the 2dseq file, suppose average=1
voxel_size = [0.167 0.167 1]; % [x y z] pixel spacing in mm
f1 = fopen(file,'r');
A=fread(f1,'int16');
A=reshape(A,dimension(1),dimension(2),dimension(3),num_dif);

% slope is a dimension(3) x num_dif matrix that defines the scaling parameters
% input the slope here
slope=[ 71.6181303266408 71.6181303266408 71.6181303266408 71.6181303266408 ...
        71.6181303266408 71.6181303266408 71.6181303266408 71.6181303266408 ...
        71.6181303266408 71.6181303266408 71.6181303266408];

slope = slope./min(min(slope));
slope = reshape(slope',dimension(3),num_dif);
for i = 1:dimension(3)
    for j = 1:num_dif
        A(:,:,i,j) = A(:,:,i,j)/slope(i,j);
    end
end

% image0 ~ imageN, where N is the number of diffusion gradient encoding
image0 = reshape(A(:,:,:,1),prod(dimension),1);
image1 = reshape(A(:,:,:,2),prod(dimension),1);
image2 = reshape(A(:,:,:,3),prod(dimension),1);
image3 = reshape(A(:,:,:,4),prod(dimension),1);
image4 = reshape(A(:,:,:,5),prod(dimension),1);
image5 = reshape(A(:,:,:,6),prod(dimension),1);
image6 = reshape(A(:,:,:,7),prod(dimension),1);
image7 = reshape(A(:,:,:,8),prod(dimension),1);
image8 = reshape(A(:,:,:,9),prod(dimension),1);
image9 = reshape(A(:,:,:,10),prod(dimension),1);
image10 = reshape(A(:,:,:,11),prod(dimension),1);

% if averaging on the repeat images is needed, run
% repeat = 5;
% dimension(3) = dimension(3)/repeat;
% image0 = reshape(sum(reshape(image0,dimension(1),dimension(2),dimension(3),repeat),4),prod(dimension),1)/repeat;
% image1 = reshape(sum(reshape(image1,dimension(1),dimension(2),dimension(3),repeat),4),prod(dimension),1)/repeat;
% image2 = reshape(sum(reshape(image2,dimension(1),dimension(2),dimension(3),repeat),4),prod(dimension),1)/repeat;
% image3 = reshape(sum(reshape(image3,dimension(1),dimension(2),dimension(3),repeat),4),prod(dimension),1)/repeat;
% image4 = reshape(sum(reshape(image4,dimension(1),dimension(2),dimension(3),repeat),4),prod(dimension),1)/repeat;
% image5 = reshape(sum(reshape(image5,dimension(1),dimension(2),dimension(3),repeat),4),prod(dimension),1)/repeat;
% image6 = reshape(sum(reshape(image6,dimension(1),dimension(2),dimension(3),repeat),4),prod(dimension),1)/repeat;

% the b-table should be arranged in the format of [b-value1 bx by bz; b-value2 bx by bz;...;b-valueN bx by bz]
b_table=[0 0 0 0 ...
1000 0.337968743631421 0.337968743631421 0
1000 0 0.337968743631421 0.337968743631421
1000 0.337968743631421 0 0.337968743631421
1000 -0.337968743631421 0.337968743631421 0
1000 0 -0.337968743631421 0.337968743631421
1000 0.337968743631421 0 -0.337968743631421
]';% [b_value bx by bz; b_value bx by bz; ...]
b_table(4,:) = -b_table(4,:);    % note that the z-dimension should be flipped for images from bruker
fclose(f1);
clear A;
clear ans;
clear f1;
save '1.src' '-v4'
%}
