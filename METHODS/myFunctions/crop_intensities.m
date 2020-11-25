function img = crop_intensities(img)
    img(img>255) = 255; 
    img(img<0) = 0;