%% Initialization
% read image
im = imread('images/chessboard03.png');
% init window for derivative in x-direction
dx = [-1 0 1;
    -1 0 1;
    -1 0 1];
% init window for derivative in y-direction
dy = dx';
% take derivatives of the image
Ix = conv2(double(im), dx, 'same');
Iy = conv2(double(im), dy, 'same');
% smoothing
gaussian = fspecial('gaussian', 9, 2);
Ix2 = conv2(Ix.*Ix, gaussian, 'same');
Iy2 = conv2(Iy.*Iy, gaussian, 'same');
Ixy = conv2(Ix.*Iy, gaussian, 'same');

% declare corner response matrices
E = zeros(size(im,1), size(im,2));
R = zeros(size(im,1), size(im,2));
% displacement from one pixel to its 8 neighbours
wx = [1 0 -1 0 1 -1 1 -1];
wy = [0 1 0 -1 1 -1 -1 1];
% constant for calculating approximation matrix R
k = 0.04;

%% Part #1,#2 - Computing matrix E with the smallest eigenvalue of M and matrix R (approximation)
tic
% For every pixel
for i = 1:size(im, 1)
    for j = 1:size(im, 2)
        m11 = 0; m12 = 0; m22 = 0;
        % calculate matrix M - local auto correlation function
        for t = 1:length(wx)
            if i+wx(t) > 0 && i+wx(t) <= size(im, 1) && j+wy(t) > 0 && j+wy(t) <= size(im, 2)
                m11 = m11 + Ix2(i+wx(t), j+wy(t));
                m12 = m12 + Ixy(i+wx(t), j+wy(t));
                m22 = m22 + Iy2(i+wx(t), j+wy(t));
            end;
        end;
        M = [m11, m12; m12, m22];
        %lambda = eig(M); % find eigenvalues of M
        %E(i,j) = min(lambda); % get minimum eigenvalue of M
        % Calculate the approximation
        E(i,j) = m11*m22-m12*m12-k*((m11+m22)^2);
        % R(i,j) = det(M)-k*((trace(M))^2);
    end;
end;
toc
figure, imshow(mat2gray(E)), title('Part 1-2: Matrix E');
Ec = E; % Copy of E, needed in other parts
%% Part #3 - Selecting 81 the most salient points (without non-maximal suppression)
n = 0; % feature counter
% declare struct for storing pixel coordinates (x, y) and values (v)
features = repmat(struct, 1, size(E,1)*size(E,2));
% convert matrix E into struct keeping pixel coordinates
for i = 1:size(im, 1)
    for j = 1:size(im, 2)
        n = n + 1;
        features(n).x = i;
        features(n).y = j;
        features(n).v = E(i,j);
    end;
end;
% sort features in descending order by value (field v)
sorted = sortStruct(features);
% display first 81 salient points
figure, imagesc(Ec), colormap(gray), hold on, axis equal, title('Part 3: Selecting 81 the most salient points (without non-maximal suppression)');
for i = 1:81
    if sorted(i).v > 0
        plot(sorted(i).y, sorted(i).x, 'gd', 'MarkerSize',5, 'MarkerFaceColor','r');
    end;
end;
hold off;
%% Part #4 Non-maximal suppression
ws = 11; % window size
for i = 1:size(im,1)-ws+1
    for j = 1:size(im,2)-ws+1
        localMax = -1; % init localMax
        x = 0; y = 0;  % and coordinates
        % find the maximum and its coordinates in current window
        for k = i:i+ws-1
            for t = j:j+ws-1
                if E(k, t) > localMax
                    localMax = E(k, t);
                    x = k; y = t;
                end;
            end;
        end;
        % set all values to zero except the maximum
        for k = i:i+ws-1
            for t = j:j+ws-1
                if ~(k == x && t == y)
                    E(k,t) = 0;
                end;
            end;
        end;
    end;
end;

clear features sorted;
% convert the matrix E to features structure
n = 0;
features = repmat(struct, 1, size(E,1)*size(E,2));
for i = 1:size(im,1)
    for j = 1:size(im,2)
        n = n + 1;
        features(n).x = i;
        features(n).y = j;
        features(n).v = E(i,j);
    end;
end;
sorted = sortStruct(features);
% display first 81 feature points
figure, imagesc(im), colormap(gray), hold on, axis equal, title('Part 4: Non-maximal suppression');
co = repmat(struct, 1, 81); % save 81 found corners: this is used in subpixel accuracy
for i = 1:81
    plot(sorted(i).y, sorted(i).x, 'gd', 'MarkerSize', 5, 'MarkerFaceColor','k');
    co(i).x = sorted(i).x; co(i).y = sorted(i).y; co(i).v = sorted(i).v;
end;
hold off;
%% Part #5 Subpixel accuracy
% Try to fit a*x^2+b*y^2+c*x*y+d*x+e*y = f
figure, imagesc(im), colormap(gray), hold on, axis equal, title('Part 5: Subpixel accuracy');
%for every corner pixel (total: 81) solve Ax=b.
for i = 1:length(co)
    A = zeros(9, 5); % init matrix A
    b = zeros(9, 1); % init vector b
    % pixels in the current window (3x3)
    tx = zeros(1, 9);
    ty = zeros(1, 9);
    tx(1) = co(i).x; ty(1) = co(i).y; b(1) = co(i).v;
    % generate matrix A:
    % add i-th pixel
    % a,            b,              c,                  d,            e
    A(1,1)=tx(1)^2; A(1,2)=ty(1)^2; A(1,3)=tx(1)*ty(1); A(1,4)=tx(1); A(1,5)=ty(1);
    % add neighbours of i-th pixel
    for j = 1:8
        tx(j+1) = tx(1)+wx(j);
        ty(j+1) = ty(1)+wy(j);
        b(j+1) = Ec(tx(j+1), ty(j+1));
        A(j+1,1)=tx(j+1)^2; A(j+1,2)=ty(j+1)^2; A(j+1,3)=tx(j+1)*ty(j+1); A(j+1,4)=tx(j+1); A(j+1,5)=ty(j+1);
    end;
    % find parameters p with least mean squares method
    %p = (A'*A)^-1*A'*b
    p = A\b;
    % get minimum coordinates
    mn = -(([2*p(1) p(3); p(3) p(2)])^-1)*([p(4); p(5)]);
    % normalize
    mn = mod(mn, 0.5);
    % plot the result
    plot(ty(1)+mn(2), tx(1)+mn(1), 'gd', 'MarkerSize', 5, 'MarkerFaceColor','k');
end;
hold off;


