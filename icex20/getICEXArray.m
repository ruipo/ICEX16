function[p,elev,az] = getICEXArray(orientation)

if strcmp(orientation,'v')

    %% Element Locations
    N = 32;
    % 0.75m spacing in the middle
    N2 = 22;
    z2 = zeros(N2,1);
    d = 0.75;
    for n = 0:N2-1
        z2(n+1) = -(n-(N2-1)/2)*d; 
    end

    % 1.5m spacing on top
    N1 = 5;
    z1 = zeros(N1,1);
    d = 1.5;
    z1(end) = z2(1) + d;

    for n = 1:N1-1
        z1(end-n) = z1(end-(n-1)) + d;
    end

    % 1.5m spacing on bottom
    N3 = 5;
    z3 = zeros(N3,1);
    d = 1.5;
    z3(1) = z2(end) - d;

    for n = 2:N3
        z3(n) = z3(n-1) - d;
    end

    z = [z1; z2; z3];
    p = [zeros(1,N) ; zeros(1,N) ; z'];
    p = p';
    
    elev = -90:1:90;
    az = 0;

elseif strcmp(orientation,'h')
    % Element Locations
    N = 32;
    % 0.75m spacing in the middle
    N2 = 22;
    z2 = zeros(N2,1);
    d = 0.75;
    for n = 0:N2-1
        z2(n+1) = -(n-(N2-1)/2)*d; 
    end

    % 1.5m spacing on top
    N1 = 5;
    z1 = zeros(N1,1);
    d = 1.5;
    z1(end) = z2(1) + d;

    for n = 1:N1-1
        z1(end-n) = z1(end-(n-1)) + d;
    end

    % 1.5m spacing on bottom
    N3 = 5;
    z3 = zeros(N3,1);
    d = 1.5;
    z3(1) = z2(end) - d;

    for n = 2:N3
        z3(n) = z3(n-1) - d;
    end

    z = [z1; z2; z3];
    p = [z'; zeros(1,N) ; zeros(1,N)]; %horizontal array
    p = p';
    
    elev = 0;
    az = 0:1:180;
    
else
    disp("Invalid orientation parameter; enter 'v' for vertical or 'h' for horizontal.")
end