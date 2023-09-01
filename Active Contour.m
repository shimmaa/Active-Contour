clear;
orgImage = imread('square.jpg');
if size(orgImage, 3) == 3
    gryImage = rgb2gray(orgImage);
else
    gryImage=orgImage;
end
%Get input from user if 1 let user determine initial contour points else initialize it as a circle 
prompt = 'Please, Enter 1 to choose initial contour or 2 to automaticly initalize it as a circle : ';
in = input(prompt);
prompt = 'Please, Enter 1 to do beta relaxation : ';
br = input(prompt);
figure
imshow (gryImage);
if in==1
    but = 1;
    i=1;
    while but == 1
        [x, y, but] = ginput(1);
        if(but==1)
            point(2,i) = x;
            point(1,i) = y;
            if  point(1,i)>size(gryImage,1)-1
                point(1,i)=size(gryImage,1)-1; 
            end
            if  point(2,i)>size(gryImage,2)-1
                point(2,i)=size(gryImage,2)-1; 
            end    
            if  point(1,i)<2
                point(1,i)=2; 
            end
            if  point(2,i)<2
                point(2,i)=2; 
            end 
            hold on
            plot(point(2,:), point(1,:), 'g-','linewidth',2);
            plot(point(2,:), point(1,:), 'y+','linewidth',2);
            i=i+1;
        else
            n=i-1;
        end
    end
    plot([point(2,end) point(2,1)], [point(1,end) point(1,1)], 'g-', 'LineWidth', 2);
elseif in==2
    % Plot the snake as a circle around the object in the image
    % Circle radius
    r = .45*min(size(gryImage,1),size(gryImage,2));
    % Cicle center [x y]
    c = [floor(size(gryImage,1)/2) floor(size(gryImage,2)/2)];
    % Number of points in the snake
    n = 80;
    % Calculate snake points in a circle
    point(1, :) = c(1)*ones(1,n)+floor(r*sin((1:n)*2*pi/n)+0.5);
    point(2, :) = c(2)*ones(1,n)+floor(r*cos((1:n)*2*pi/n)+0.5);
    hold on
    plot(point(2,:), point(1,:), 'g-','linewidth',2);
    plot(point(2,:), point(1,:), 'y+','linewidth',2);
    plot([point(2,end) point(2,1)], [point(1,end) point(1,1)], 'g-', 'LineWidth', 2);
end
if ~isempty(in)
    if in(1)==1 || in(1)==2 
        alpha = 0.05; beta = 1; gamma = 5; maxIt =100;
        % Set the curvature threshold
        cThreshold = 0.1;
        % Set the image energy threshold
        imgEnrgT = 1;
        % Round indices of snake points
        point(1:2,:) = round(point(1:2,:));
        prevpoint(1:2,:) = circshift(point(1:2,:),[0 1]);
        nextpoint(1:2,:) = circshift(point(1:2,:),[0 -1]);
        point(3,:) = alpha;
        point(4,:) = beta;
        point(5,:) = gamma;
        %Computes averge distance of the snake
        sum=0;
        for i = 1:n
            sum = sum + norm(nextpoint(1:2, i)- point(1:2, i));
        end
        avgDist = sum/n;
        % Get the image energy
%         % Build Gauss filter
%         sigma = 3;
%         fSize = [ceil(3*sigma) ceil(3*sigma)];
%         gf = fspecial('gaussian', fSize, sigma);
%         % Apply Gauss filter to image
%         gaussImg = double(imfilter(gryImage, gf, 'replicate'));
%         % Find gradients using a Sobel filter
%         fy = fspecial('sobel');
%         fx = fy';
%         Iy = imfilter(gaussImg, fy, 'replicate');
%         Ix = imfilter(gaussImg, fx, 'replicate');
%         enrgImg = sqrt(Ix.^2 + Iy.^2);
        enrgImg = edge(gryImage,'Canny');
        %imshow(enrgImg);      
        % Main loop that evolves the snake
        it=1;
        flag=true;
        while flag==true
            % Counter for the number of points moved
            pointsMoved = 0;
            for i = 1:n
                % Extract pixels in the neighborhood of the current snake point
                neighborhood = enrgImg((point(1, i)-1):(point(1, i)+1),(point(2, i)-1):(point(2, i)+1));
                % Normalise the image energy 
                enrgMin = min(min(neighborhood));
                enrgMax = max(max(neighborhood));
                if (enrgMax - enrgMin) < 5
                    enrgMin = enrgMax - 5;
                end
                    normenrg = (enrgMin - neighborhood) / (enrgMax - enrgMin);
               for wm=-1:1:1
                   for wn=-1:1:1 
                        % Current position
                        pos = point(1:2,i) + [wm wn]';
                        % Calculate the continuity term
                        Econt(wm+2,wn+2) = abs(norm(nextpoint(1:2, i)- pos) - avgDist);
                        %Econt(wm+2,wn+2) = norm(nextpoint(1:2, i)- pos);
                        % Calculate the curvature term
                        Ecurv(wm+2,wn+2) = norm(prevpoint(1:2, i) - 2*pos + nextpoint(1:2,i))^2;
                   end
               end
                % Normalize the continuity and curvature terms to lie in the range [0,1]
                Ecurv = Ecurv / max(max(Ecurv));
                Econt = Econt / max(max(Econt));
                Eimg = normenrg;
                % Sum the energy terms
                En = point(3,i)*Econt + point(4,i)*Ecurv + point(5,i)*Eimg;
                % Find the location of minimum energy in the neighborhood
                [r, c] = find(En == min(min(En)));
                point(6,i) = neighborhood(r(1),c(1));
                r=r(1)-2;
                c=c(1)-2;
                if r~=0 || c~=0
                    if point(1,i)>2 && point(1,i)<size(gryImage,1)-1 && point(2,i)>2 && point(2,i)<size(gryImage,2)-1
                        point(1:2,i) = point(1:2, i) + [r c]';
                        % Increment counter
                        pointsMoved = pointsMoved + 1;  
                    end
                end
            end
            if (it == maxIt || pointsMoved < 3)
                flag = false;
            else
                prevpoint(1:2,:) = circshift(point(1:2,:),[0 1]);
                nextpoint(1:2,:) = circshift(point(1:2,:),[0 -1]);
                if br==1
                    % Iterate through all snake points to find curvatures
                    for i = 1:n
                    % Estimate the curvature at each point
                    ui = point(1:2,i) - prevpoint(1:2,i);
                    uiPlus = nextpoint(1:2,i) - point(1:2,i);
                    cr(i) = norm( ui/norm(ui) - uiPlus/norm(uiPlus) )^2;
                    end
                    % Iterate through all snake points to find where to relax beta
                    pri = circshift(1:n,[0 1]);
                    nxi = circshift(1:n,[0 -1]);
                    for i = 1:n
                        if cr(i)>cr(pri(i)) && cr(i)>cr(nxi(i)) && cr(i)>cThreshold && point(6,i)>=imgEnrgT && point(4,i)~=0
                            point(4,i)=0;
                        disp(['Relaxed beta for point nr. ', num2str(i)]);
                        end
                    end
                end
            end
            it=it+1;
        end
        %Draw final contour
        hold on
        plot(point(2,:), point(1,:), 'r-','linewidth',2);
        plot(point(2,:), point(1,:), 'b+','linewidth',2);
        plot([point(2,end) point(2,1)], [point(1,end) point(1,1)], 'r-', 'LineWidth', 2);
    end
end