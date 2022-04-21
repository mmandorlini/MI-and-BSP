function FBP_comparison_D(image, reference)

n_pixel = size(image,1); % image dimension in pixels
pixel_size = 1; % pixel dimension in mm
fs = 1/pixel_size; % spatial sampling frequency
delta_f = fs/n_pixel;
f_axis = [0:delta_f:fs-delta_f]; 

N = [180]; % number of projections angles
D = [0.2:0.2:1]; % fraction of frequencies below the nyquist which we want to pass
filter_name = ["none", "ram-lak", "shepp-logan", "cosine", "hamming", "hann"]; % type of filter


% iteration using different numbers of projection angles
for n = 1:length(N) 
    
    MSE = []; RMSE = []; 

    
    % iteration using different fraction of frequencies
    for di = 1:length(D)
        
        figure, sgtitle(['d = ', num2str(D(di))]);
        
        % iteration using different filters   
        for f = 1:length(filter_name) 
            
% [[[SINOGRAM]]]             
            theta = linspace(0,180,N(n)); % angles vector
            num_angles = length(theta); % number of angles
            parallel_proj = size(image,1); % number of parallel projections lines
            
            sinogram = zeros(parallel_proj,num_angles); % matrix to store parallel projections
            
            % creation of the sinogram
            for i = 1:num_angles
                imgrotation = imrotate(image, 90-theta(i), 'bilinear', 'crop');  % rotation of the image of an angle theta(i)
                sinogram(:,i) = (sum(imgrotation))';  % line integration on the rotated image columns
            end
            
% [[[FILTERING]]] 
            if filter_name(f) == "none" % no filtering case
            else              
                filter = filter_name(f);      
                d = D(di);
                len = n_pixel/2;
                H = designFilter(filter, len, d); % definition of the high pass filter
                
                sinogram = fft(sinogram);  % FFT to obtain sinogram in frequency domain
                
                for i=1:num_angles
                    sinogram(:,i) = sinogram(:,i).*H; % actual filtering of each projection in frequency domain
                end
                
                sinogram = ifft(sinogram,'symmetric'); % inverse FFT to obtain sinogram in spatial domain 
            end
            
            % filtered sinogram visualization
            subplot(4,length(filter_name),f+length(filter_name)), imagesc(sinogram);

% [[[RECONSTRUCTION]]]
            theta = (pi/180)*theta; % convert theta vector from degrees to radians
            im_rec = zeros(parallel_proj); % matrix to stora the reconstructed image

            origin = floor(parallel_proj/2) + 1;  % central point of the projection plane
            [x,y] = meshgrid(ceil(-parallel_proj/2):ceil(parallel_proj/2-1)); % coordinates of the plane

            subplot(4,length(filter_name),f),
            
            % actual image reconstruction
            for i = 1:num_angles
                t = round(origin + x*sin(theta(i)) + y*cos(theta(i))); % definition of the new polar coordinates 
                indices = find((t > 0) & (t <= parallel_proj)); % selecting the new coordinates inside the image bounds
                t = t(indices);

                contrib = sinogram(:,i); % contribution of each projection on the final image

                % summation
                im_rec(indices) = im_rec(indices) + contrib(t)./num_angles; % actual reconstruction of the image, using only values inside the image bounds

                % visualization
                imshow(im_rec,[]);
                title([filter_name(f),' filter']);
            end
            
            if filter_name == "none" % no filtering case
                im_rec = rescale(im_rec); 
            end
            
            % visualization
            subplot(4,length(filter_name),f+2*length(filter_name)), histogram(im_rec);
            Fourier_im_rec = fftshift(fft2(im_rec));
            subplot(4,length(filter_name),f+3*length(filter_name)), mesh(abs(Fourier_im_rec));
            
% [[[PERFORMANCE EVALUATION]]]
% with respect to the original image without noise
            MSE(f) = immse(im_rec,reference); % mean square error
            RMSE(f) = sqrt(MSE(f)); % root mean square error
            PSNR(f) = psnr(im_rec, reference); % peak signal to noise ratio
            SSIM(f) = ssim(im_rec, reference); % structural similarity

        end
        
        % display evaluation metrics
        disp('MSE:')
        fprintf([num2str(MSE(1),2), num2str(MSE(1:6),'    %.4f'),'\n\n'])
        disp('RMSE:')
        disp(RMSE)
        disp('PSNR:')
        disp(PSNR)
        disp('SSIM:')
        disp(SSIM)
    end 
end 
end