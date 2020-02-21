function mov1 = mov_resize(mov, siz)
% mov: 4D uint8
% siz: [height, width]

mov1 = zeros(siz(1), siz(2), size(mov,3), size(mov,4), 'like', mov);

for iF = 1:size(mov,4)
    mov1(:,:,:,iF) = imresize(mov(:,:,:,iF), siz, 'Colormap', 'original');
end 
end %func