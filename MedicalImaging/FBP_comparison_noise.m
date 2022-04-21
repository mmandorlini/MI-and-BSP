function FBP_comparison_noise(image, reference)

n_pixel = size(image,1);
im_dim = n_pixel; % image dimension
pixel_size = 0.001*10^3; % pixel dimension from m to mm
fs = 1/pixel_size; % spacial sampling frequency
delta_f = fs/n_pixel;
f_axis = [0:delta_f:fs-delta_f];

N = [10,50,100,180]; % number of projections angles
D = [1]; % fraction of frequencies below the nyquist which we want to pass
filter_name = ["none", "ram-lak", "shepp-logan", "cosine", "hamming", "hann"]; %different types of filters


%iteration using different numbers of projection angles____________________
for n = 1:length(N) 
    
    MSE = []; RMSE = []; 
    figure, sgtitle(['N=', num2str(N(n))]);
    %iteration using different fraction of frequencies_____________________
    for di = 1:length(D)
        
        %iteration using different filters_________________________________   
        for f = 1:length(filter_name) 
            
            % [[[SINOGRAM]]]  
            theta = linspace(0,180,N(n)); % angles vector
            num_angles=length(theta); % number of angles

            parallel_proj=size(image,1); % number of parallel projections lines
            sinogram = zeros(parallel_proj,num_angles); % matrix to store parallel projections
            
            % creation of the sinogram
            for i = 1:num_angles
                imgrotation = imrotate(image, 90-theta(i), 'bilinear', 'crop');
                sinogram(:,i) = (sum(imgrotation))';  
            end
            
            % [[[FILTERING]]]
            if filter_name(f) == "none" % no filtering
                sinogram;
            else              
            % sinogram filtering in frequency domain
                filter = filter_name(f);      
                d=D(di);
                len = n_pixel/2;
                H = designFilter(filter, len, d); 
                sinogram = fft(sinogram); 
                
                % actual filtering of each projection in frequency domain
                for i=1:num_angles
                    sinogram(:,i) = sinogram(:,i).*H; 
                end
                
                sinogram = ifft(sinogram,'symmetric');
            end
        
            subplot(4,length(filter_name),f+length(filter_name)), imagesc(sinogram);
    
            % [[[RECONSTRUCTION]]]
            theta = (pi/180)*theta; % convert theta to radians
            im_rec = zeros(parallel_proj); % storage of image reconstructions

            midindex = floor(parallel_proj/2) + 1; % middle index of the projection lines
            [xCoords,yCoords] = meshgrid(ceil(-parallel_proj/2):ceil(parallel_proj/2-1)); % set up the coords of the image

            subplot(4,length(filter_name),f),
            
            % algorithm for image reconstruction
            for i = 1:num_angles
                % selecting the coordinates inside the bounds
                Coords = round(midindex + xCoords*sin(theta(i)) + yCoords*cos(theta(i)));
                indices = find((Coords > 0) & (Coords <= parallel_proj)); 
                Coords = Coords(indices);

                contrib = sinogram(:,i);

                % summation
                im_rec(indices) = im_rec(indices) + contrib(Coords)./num_angles;

                % visualization
                imagesc(im_rec), colormap(gray);
                drawnow, title([filter_name(f),' filter']);
            end
            
            if filter_name == "none" % in case of no filtering
                im_rec = rescale(im_rec);
            end
            
            subplot(4,length(filter_name),f+2*length(filter_name)), histogram(im_rec);
            Fourier_im_rec = fftshift(fft2(im_rec));
            subplot(4,length(filter_name),f+3*length(filter_name)), mesh(abs(Fourier_im_rec));

            % [[[PERFORMANCE EVALUATION]]]
            MSE(f) = immse(im_rec,reference); % mean square error
            RMSE(f) = sqrt(MSE(f)); % root mean square error
            PSNR(f) = psnr(im_rec, reference); % peak signal to noise ratio
            SSIM(f) = ssim(im_rec, reference); % structural similarity
        
            % fprintf([filter_name, ' filter parameters:\n'])
            % fprintf(['MSE = ',num2str(MSE(f)),'  RMSE = ',num2str(RMSE(f)),'  PSNR = ',num2str(PSNR(f)),'  SSIM = ',num2str(SSIM(f)),])

        end %_________________________end iteration using different filters
        
        % display evaluation metrics
        disp('MSE:')
        disp(MSE)
        disp('RMSE:')
        disp(RMSE)
        disp('PSNR:')
        disp(PSNR)
        disp('SSIM:')
        disp(SSIM)
    end %___________________end iteration using different fraction of freqs

end %__________________end iteration using different numbers of proj angles
end