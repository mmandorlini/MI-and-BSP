function [im_rec, sinogram] = Image_Reconstruction_noise(image, N, filter_name, d)

reference = phantom('Modified Shepp-Logan'); % set the reference image ( original phantom without noise superimposition)

figure,
% [[[SINOGRAM]]]
    theta = linspace(0,180,N); % angles vector
    n_pixel = size(image,1); % image size in pixels 
    parallel_proj = size(image,1); % number of parallel projections
    sinogram = zeros(parallel_proj,N); % matrix to store parallel projections

    for i = 1:N
        imgrotation = imrotate(image, 90-theta(i), 'bilinear', 'crop');  % rotation of the image of an angle theta(i)
        sinogram(:,i) = (sum(imgrotation))';  % line integration on the rotated image columns
    end
    
    figure,    
    subplot(323), imshow(sinogram,[]);
        
% [[[FILTERING]]]
   if filter_name == "none" % no filtering case
   else 
        filter = filter_name;    
        len = n_pixel/2;   
        H = designFilter(filter, len, d);  % definition of the high pass filter 
    
        sinogram = fft(sinogram); % FFT to obtain sinogram in frequency domain
    
        for i=1:N
            sinogram(:,i) = sinogram(:,i).*H; % actual filtering of each projection in frequency domain
        end
    
        sinogram = ifft(sinogram,'symmetric'); % inverse FFT to obtain sinogram in spatial domain 
    end
    
% [[[RECONSTRUCTION]]]
    theta = (pi/180)*theta; % convert theta vector from degrees to radians
    im_rec = zeros(parallel_proj); % storage of image reconstructions

     origin = floor(parallel_proj/2) + 1; % central point of the projection plane
    [x,y] = meshgrid(ceil(-parallel_proj/2):ceil(parallel_proj/2-1)); % coordinates of the plane

    for i = 1:N

        % selecting the coordinates inside the bounds
        t = round(origin + x*sin(theta(i)) + y*cos(theta(i)));
        indices = find((t > 0) & (t <= parallel_proj)); 
        t = t(indices);
        
        contrib = sinogram(:,i);
    
        % summation
        im_rec(indices) = im_rec(indices) + contrib(t)./N;
    
    end

% [[[PERFORMANCE EVALUATION]]]
% with respect to the original image without noise
    MSE = immse(im_rec, reference);
    RMSE = sqrt(MSE);
    PSNR = psnr(im_rec, reference);
    SSIM = ssim(im_rec, reference);

    % display evaluation metrics
    disp(['MSE: ', num2str(MSE)])
    disp(['RMSE: ', num2str(RMSE)])
    disp(['PSNR: ', num2str(PSNR)])
    disp(['SSIM: ', num2str(SSIM)])

% [[[VISUALIZATION]]]
    subplot(331), imshow(reference,[]);
    subplot(332), imshow(image,[]);
    subplot(332), imshow(im_rec,[]);
    subplot(333), histogram(reference);
    subplot(334), histogram(image);
    subplot(335), histogram(im_rec);
    Fourier = fftshift(fft2(reference));
    subplot(336), mesh(abs(Fourier));
    Fourier = fftshift(fft2(image));
    subplot(337), mesh(abs(Fourier));
    Fourier = fftshift(fft2(im_rec));
    subplot(338), mesh(abs(Fourier));
end