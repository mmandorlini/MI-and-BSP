classdef Functions
% Collection of functions necessary for code : 
properties
end
methods (Static)
    
function [MSE,RMSE,PSNR,SSIM] = FBP_comparison(image, reference)

% create figures for the comparison of different FBP algorithms in case of
% reconstruction of an image and print some reconstruction quality indices

% reconstruction are obtained varing:
% N = number of angles of projection
% filter_name = type of filter applied to the projection

% reconstructions are comparerd in terms of:
% MSE = mean square error
% RMSE = root mean square error
% PSNR = peak signal to noise ratio
% SSIM = image structural similarity

% Input: 
%       image = image to be reconstructed
%       reference = image of reference for the recunstruction quality
%       evaluation


 

n_pixel = size(image,1); % image dimension in pixels
pixel_size = 1; % pixel dimension in mm
fs = 1/pixel_size; % spatial sampling frequency
delta_f = fs/n_pixel;
f_axis = [0:delta_f:fs-delta_f]; 

N = [10,50,100,180]; % number of projections angles
D = [1]; % fraction of frequencies below the nyquist which we want to pass
filter_name = ["none", "ram-lak", "shepp-logan", "cosine", "hamming", "hann"]; % type of filter


% iteration using different numbers of projection angles
for n = 1:length(N) 

    figure, sgtitle(['N = ', num2str(N(n))]);
    
    % iteration using different fraction of frequencies
    for di = 1:length(D)
         
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
                H = Functions.designFilter(filter, len, d); % definition of the high pass filter
                
                sinogram = fft(sinogram);  % FFT to obtain sinogram in frequency domain
                
                for i=1:num_angles
                    sinogram(:,i) = sinogram(:,i).*H; % actual filtering of each projection in frequency domain
                end
                
                sinogram = ifft(sinogram,'symmetric'); % inverse FFT to obtain sinogram in spatial domain 
            end
            
            % filtered sinogram visualization
            subplot(4,length(filter_name),f+length(filter_name)), imagesc(sinogram);
    
% [[[RECONSTRUCTION]]]
            theta_deg = theta;
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
            
            % filtered sinogram visualization
            subplot(4,length(filter_name),f+length(filter_name)), imagesc(theta_deg,t-origin,sinogram), xlabel('theta'),ylabel('t');
            
            if filter_name == "none" % no filtering case
                im_rec = rescale(im_rec); 
            end
            
            % visualization
            subplot(4,length(filter_name),f+2*length(filter_name)), histogram(im_rec);
            Fourier_im_rec = fftshift(fft2(im_rec));
            subplot(4,length(filter_name),f+3*length(filter_name)), mesh(abs(Fourier_im_rec));

% [[[PERFORMANCE EVALUATION]]] 
% with respect to the original image without noise
            MSE(n,f) = immse(im_rec,reference); % mean square error
            RMSE(n,f) = sqrt(MSE(f)); % root mean square error
            PSNR(n,f) = psnr(im_rec, reference); % peak signal to noise ratio
            SSIM(n,f) = ssim(im_rec, reference); % structural similarity
        
        end 
        
        % display evaluation indices        
        disp('MSE:')
        fprintf([num2str(MSE(n,1),2), num2str(MSE(n,1:6),'    %.4f'),'\n\n'])
        disp('RMSE:')
        disp(RMSE(n,:))
        disp('PSNR:')
        disp(PSNR(n,:))
        disp('SSIM:')
        disp(SSIM(n,:))
    end 
end 
end


function [MSE,RMSE,PSNR,SSIM] = FBP_comparison_D(image, reference)

% create figures for the comparison of different FBP algorithms in case of
% reconstruction of an image and print some reconstruction quality indices

% reconstruction are obtained varing:
% filter_name = type of filter applied to the projection
% D = fraction of frequencies below the Nyquist which we want to pass

% reconstructions are comparerd in terms of:
% MSE = mean square error
% RMSE = root mean square error
% PSNR = peak signal to noise ratio
% SSIM = image structural similarity

% Input: 
%       image = image to be reconstructed
%       reference = image of reference for the recunstruction quality
%       evaluation


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
                H = Functions.designFilter(filter, len, d); % definition of the high pass filter
                
                sinogram = fft(sinogram);  % FFT to obtain sinogram in frequency domain
                
                for i=1:num_angles
                    sinogram(:,i) = sinogram(:,i).*H; % actual filtering of each projection in frequency domain
                end
                
                sinogram = ifft(sinogram,'symmetric'); % inverse FFT to obtain sinogram in spatial domain 
            end
            

% [[[RECONSTRUCTION]]]
           theta_deg = theta;
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
            
            % filtered sinogram visualization
            subplot(4,length(filter_name),f+length(filter_name)), imagesc(theta_deg,t-origin,sinogram), xlabel('theta'),ylabel('t');
            
            if filter_name == "none" % no filtering case
                im_rec = rescale(im_rec); 
            end
            
            % visualization
            subplot(4,length(filter_name),f+2*length(filter_name)), histogram(im_rec);
            Fourier_im_rec = fftshift(fft2(im_rec));
            subplot(4,length(filter_name),f+3*length(filter_name)), mesh(abs(Fourier_im_rec));
            
% [[[PERFORMANCE EVALUATION]]]
% with respect to the original image without noise
            MSE(di,f) = immse(im_rec,reference); % mean square error
            RMSE(di,f) = sqrt(MSE(f)); % root mean square error
            PSNR(di,f) = psnr(im_rec, reference); % peak signal to noise ratio
            SSIM(di,f) = ssim(im_rec, reference); % structural similarity

        end
        
        % display evaluation metrics
        disp('MSE:')
        fprintf([num2str(MSE(di,1),2), num2str(MSE(di,1:6),'    %.4f'),'\n\n'])
        disp('RMSE:')
        disp(RMSE(di,:))
        disp('PSNR:')
        disp(PSNR(di,:))
        disp('SSIM:')
        disp(SSIM(di,:))
    end 
end 
end


function filt = designFilter(filter, len, d)
% Returns the Fourier Transform of the filter which will be
% used to filter the projections
%
% INPUT ARGS:   filter - either the string specifying the filter
%               len    - the length of the projections
%               d      - the fraction of frequencies below the nyquist
%                        which we want to pass
%
% OUTPUT ARGS:  filt   - Fourier transform of the filter to use on the projections


order = max(64,2^nextpow2(2*len));

if strcmpi(filter, 'none')
    filt = ones(1, order);
    return;
end

% First create a bandlimited ramp filter (Eqn. 61 Chapter 3, Kak and
% Slaney) - go up to the next highest power of 2.

n = 0:(order/2); % 'order' is always even. 
filtImpResp = zeros(1,(order/2)+1); % 'filtImpResp' is the bandlimited ramp's impulse response (values for even n are 0)
filtImpResp(1) = 1/4; % Set the DC term 
filtImpResp(2:2:end) = -1./((pi*n(2:2:end)).^2); % Set the values for odd n
filtImpResp = [filtImpResp filtImpResp(end-1:-1:2)]; 
filt = 2*real(fft(filtImpResp)); 
filt = filt(1:(order/2)+1);

w = 2*pi*(0:size(filt,2)-1)/order;   % frequency axis up to Nyquist

switch filter
    case 'ram-lak'
        % Do nothing
    case 'shepp-logan'
        % be careful not to divide by 0:
        filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));
    case 'cosine'
        filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));
    case 'hamming'
        filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));
    case 'hann'
        filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./d)) / 2;
    otherwise
        error(message('images:iradon:invalidFilter'))
end

filt(w>pi*d) = 0;                      % Crop the frequency response
filt = [filt' ; filt(end-1:-1:2)'];    % Symmetry of the filter

end


function [best_filter, best_N, best_D, best_MSE] = MSE_Optimal_Filter(image, reference)
% returns the parameter combination for the best filter in terms of MSE in
% a predefined ROI (square region from pixel (200,200) to pixel (300,300)).

% reconstruction are obtained varing:
% N = number of angles of projection
% filter_name = type of filter applied to the projection
% D = fraction of frequencies below the Nyquist which we want to pass

% reconstructions are comparerd in terms of:
% MSE = mean square error

% Input: 
%       image = image to be reconstructed
%       reference = image of reference for the recunstruction quality
%       evaluation


n_pixel = size(image,1); % image dimension in pixels
pixel_size = 1; % pixel dimension in mm
fs = 1/pixel_size; % spatial sampling frequency
delta_f = fs/n_pixel;
f_axis = [0:delta_f:fs-delta_f]; 

N = [180]; % number of projections angles
D = [0.2:0.05:1]; % fraction of frequencies below the nyquist which we want to pass
filter_name = ["none", "ram-lak", "shepp-logan", "cosine", "hamming", "hann"]; % type of filter

% iteration using different numbers of projection angles
for n = 1:length(N)
    
    MSE = [];

% iteration using different fraction of frequencies
for di = 1:length(D)
    
    % iteration using different filters 
    for f = 1:length(filter_name) 

% [[[SINOGRAM]]] 
        theta = linspace(0,180,N(n)); % angles vector
        num_angles = length(theta); % number of angles
 
        parallel_proj=size(image,1); % number of parallel projections lines
        
        sinogram = zeros(parallel_proj,num_angles); % matrix to store parallel projections
        
        % creation of the sinogram
        for i = 1:num_angles
            imgrotation = imrotate(image, 90-theta(i), 'bilinear', 'crop'); % rotation of the image of an angle theta(i)
            sinogram(:,i) = (sum(imgrotation))';  % line integration on the rotated image columns
        end

% [[[FILTERING]]] 
       if filter_name(f) == "none" %no filtering case
       else  
            % sinogram filtering in freq domain
            filter = filter_name(f);      
            d=D(di);
            len = n_pixel/2;
            H = Functions.designFilter(filter, len, d); % definition of the high pass filter
            
            sinogram = fft(sinogram);  % FFT to obtain sinogram in frequency domain
                
            for i=1:num_angles
                sinogram(:,i) = sinogram(:,i).*H; % actual filtering of each projection in frequency domain
            end
                
            sinogram = ifft(sinogram,'symmetric'); % inverse FFT to obtain sinogram in spatial domain 
        
       end
    
% [[[RECONSTRUCTION]]]
        theta = (pi/180)*theta; % convert theta vector from degrees to radians
        im_rec = zeros(parallel_proj); % matrix to stora the reconstructed image

        origin = floor(parallel_proj/2) + 1;  % central point of the projection plane
        [x,y] = meshgrid(ceil(-parallel_proj/2):ceil(parallel_proj/2-1)); % coordinates of the plane

        % actual image reconstruction
        for i = 1:num_angles
            t = round(origin + x*sin(theta(i)) + y*cos(theta(i))); % definition of the new polar coordinates 
            indices = find((t > 0) & (t <= parallel_proj)); % selecting the new coordinates inside the image bounds
            t = t(indices);
        
            contrib = sinogram(:,i); % contribution of each projection on the final image
    
            % summation
            im_rec(indices) = im_rec(indices) + contrib(t)./num_angles; %actual reconstruction of the image, using only values inside the image bounds
    
        end

% [[[PERFORMANCE EVALUATION]]]  with respect to che original image
        MSE(f,1) = immse(im_rec(200:300,200:300),reference(200:300,200:300));
        best_MSE_filter(di,n) = filter_name(MSE == min(MSE));

    end
    MSE_D(di) = min(MSE(:));
    best_MSE_D(n) = D(MSE_D == min(MSE_D));
    
end
MSE_N(n) = min(MSE_D(:));
best_MSE_N = N(MSE_N == min(MSE_N));  

end

% save parameters that guarantee the best performance in terms of MSE
best_N = best_MSE_N;
best_D = best_MSE_D(MSE_N == min(MSE_N));
best_filter = best_MSE_filter(best_D == best_MSE_D, MSE_N == min(MSE_N));
best_MSE = min(MSE_N);

% display evaluation indices 
disp(['best_filter: ', char(best_filter)])
disp(['best_N: ', num2str(best_N)])
disp(['best_D: ',num2str(best_D)])
disp(['best_MSE: ', num2str(best_MSE)])

end


function [im_rec, sinogram] = FBP(image, N, filter_name, d)
% performs the FBP reconstruction using the input parameters in case of
% input image without noise

% reconstructions are comparerd in terms of:
% MSE = mean square error
% RMSE = root mean square error
% PSNR = peak signal to noise ratio
% SSIM = image structural similarity

% Input:
%       image = image to be reconstructed
%       N = number of angles of projection
%       filter_name = type of filter
%       d = fraction of frequencies below the Nyquist which we want to pass 

% Output:
%       im_rec = reconstructed image
%       sinogram = sinogram obtained after the filtering phase


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
        H = Functions.designFilter(filter, len, d);  % definition of the filter in frequency domain
    
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
        t = round(origin + x*sin(theta(i)) + y*cos(theta(i))); % definition of the new polar coordinates 
        indices = find((t > 0) & (t <= parallel_proj)); % selecting the new coordinates inside the image bounds
        t = t(indices);
        
        contrib = sinogram(:,i);
    
        % summation
        im_rec(indices) = im_rec(indices) + contrib(t)./N;
    
    end

% [[[PERFORMANCE EVALUATION]]] 
% with respect to the original image
    MSE = immse(im_rec, image);
    RMSE = sqrt(MSE);
    PSNR = psnr(im_rec, image);
    SSIM = ssim(im_rec, image);

% [[[VISUALIZATION]]]
    subplot(321), imshow(image,[]), title('original image');
    subplot(322), imshow(im_rec,[]), title('reconstructed image');
    subplot(323), histogram(image);
    subplot(324), histogram(im_rec);
    Fourier = fftshift(fft2(image));
    subplot(325), mesh(abs(Fourier));
    Fourier = fftshift(fft2(im_rec));
    subplot(326), mesh(abs(Fourier));

    % display evaluation indices 
    disp(['MSE: ', num2str(MSE)])
    disp(['RMSE: ', num2str(RMSE)])
    disp(['PSNR: ', num2str(PSNR)])
    disp(['SSIM: ', num2str(SSIM)])

end


function [im_rec, sinogram] = FBP_noise(image, reference, N, filter_name, d)
% performs the FBP reconstruction using the input parameters in case of
% input image without noise

% reconstructions are comparerd in terms of:
% MSE = mean square error
% RMSE = root mean square error
% PSNR = peak signal to noise ratio
% SSIM = image structural similarity

% Input:
%       image = image to be reconstructed
%       reference = image of reference for the recunstruction quality
%       N = number of angles of projection
%       filter_name = type of filter
%       d = fraction of frequencies below the Nyquist which we want to pass 

% Output:
%       im_rec = reconstructed image
%       sinogram = sinogram obtained after the filtering phase


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
        H = Functions.designFilter(filter, len, d);  % definition of the high pass filter 
    
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
        t = round(origin + x*sin(theta(i)) + y*cos(theta(i))); % definition of the new polar coordinates 
        indices = find((t > 0) & (t <= parallel_proj)); % selecting the new coordinates inside the image bounds
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
    subplot(331), imshow(reference,[]), title('reference image');
    subplot(332), imshow(image,[]), title('original noisy image');
    subplot(333), imshow(im_rec,[]), title('reconstructed image');
    subplot(334), histogram(reference);
    subplot(335), histogram(image);
    subplot(336), histogram(im_rec);
    Fourier = fftshift(fft2(reference));
    subplot(337), mesh(abs(Fourier));
    Fourier = fftshift(fft2(image));
    subplot(338), mesh(abs(Fourier));
    Fourier = fftshift(fft2(im_rec));
    subplot(339), mesh(abs(Fourier));
end


end
end