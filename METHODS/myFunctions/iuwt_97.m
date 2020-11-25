function [w, s_out] = iuwt_97(im, levels)
% See "iuwt"
% This function computes 2D IUWT using 9/7 filterbank 
%
% "JPEG-2000 uses the biorthogonal 9/7 wavelet transform"
% "The 9/7 biorthogonal and 5/3 wavelets in Figure 7.15 are 
% recommended for image compression in JPEG-2000 and are often
% used in wavelet image-processing applications. 
% The 5/3 biorthogonal wavelet has p1=p2=2 vanishing moments, while the
% 9/7 wavelet has p1=p2=4 vanishing moments" Wavelet tours
%
% ASHKAN
%

% Determine output class - single if input is single, or contains many elements
if strcmp(class(im), 'single') || numel(im) > 10000000
    wclass = 'single';
    s_in = single(im);
else
    wclass = 'double';
    s_in = double(im);
end   

% Preallocate wavelet output; 3-d even if input is a vector
w = zeros([size(im) length(levels)], wclass);

% low-pass filter in the CDF 9/7 filters
lp = [.026748757411 -.016864118443 -.078223266529 .266864118443];
kernel = [lp .602949018236 fliplr(lp)];
% The above kernel found at 
% https://www.mathworks.com/matlabcentral/fileexchange/11846-cdf-9-7-wavelet-transform


% Compute transform
for ii = 1:levels(end)
    % Create convolution kernel
    h = dilate_wavelet_kernel(kernel, 2^(ii-1)-1);
    
    % Convolve and subtract to get wavelet level
%     s_out = imfilter(s_in, h' * h, padding);
    s_out = conv2(s_in,h'*h,'same'); %%%%%%%%%%%%%%% ASHKAN %%%%%%%%%%

    % Store wavelet level only if it's in LEVELS
    ind = find(levels == ii);
    if isscalar(ind)
        w(:,:,ind) = s_in - s_out;
    end
    
    % Update input for new iteration
    s_in = s_out;
end

% Remove singleton dimensions
w = squeeze(w);



function h2 = dilate_wavelet_kernel(h, spacing)
% Dilates a wavelet kernel by entering SPACING zeros between each
% coefficient of the filter kernel H.

% Check input
if ~isvector(h) && ~isscalar(spacing)
    error(['Invalid input to DILATE_WAVELET_KERNEL: ' ...
          'H must be a vector and SPACING must be a scalar']);
end

% Preallocate the expanded filter
h2 = zeros(1, numel(h) + spacing * (numel(h) - 1));
% Ensure output kernel orientation is the same
if size(h,1) > size(h,2)
    h2 = h2';
end
% Put in the coefficients
h2(1:spacing+1:end) = h;