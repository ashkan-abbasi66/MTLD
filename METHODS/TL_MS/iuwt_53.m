function [w, s_out] = iuwt_53(im, levels)
% See "iuwt"
% This function computes 2D IUWT using 5/3 filterbank 
% 
% Found at
% "The 5/3 and 9/7 wavelet filters study in a sub-bands image coding"
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

% coefficients for filter
kernel =[-1/8,1/4,3/4,1/4,-1/8];

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