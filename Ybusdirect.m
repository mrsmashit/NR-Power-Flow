clc;
clearvars;
%% Read Line Data

linedata = xlsread('zinfo_4bus.xls', 'Sheet1', 'A1:M9'); % this is now only wye-delta in zinfo file
nolines = length(linedata(:, 1));
fromnode = linedata(:, 1);
tonode = linedata(:, 2);

for loop1 = 1:3:nolines

    if linedata(loop1, 12) == 1 % means feeder
        zinfo = linedata(loop1:loop1 + 2, 3:5) + j * linedata(loop1:loop1 + 2, 6:8);
        binfo = j * linedata(loop1:loop1 + 2, 9:11);
        Ynn(:, :, fromnode(loop1), tonode(loop1)) = inv(zinfo);
        Ynn(:, :, tonode(loop1), fromnode(loop1)) = inv(zinfo);
        Ynm(:, :, fromnode(loop1), tonode(loop1)) = inv(zinfo) + 0.5 * binfo;
        Ynm(:, :, tonode(loop1), fromnode(loop1)) = inv(zinfo) + 0.5 * binfo;
    end

    %% grounded Wye- Grounded Wye connection
    if linedata(loop1, 12) == 2

        alpha = linedata(loop1, 13);
        zinfo = linedata(loop1:loop1 + 2, 3:5) + j * linedata(loop1:loop1 + 2, 6:8);
        Y1 = inv(zinfo);
        Ynn(:, :, fromnode(loop1), tonode(loop1)) = Y1 / (alpha * alpha);
        Ynn(:, :, tonode(loop1), fromnode(loop1)) = Y1;
        Ynm(:, :, fromnode(loop1), tonode(loop1)) = Y1 / alpha; % this sign needs be revised
        Ynm(:, :, tonode(loop1), fromnode(loop1)) = Y1 / alpha; % this sign needs be revised
    end

    %% Del-Del connection
    if linedata(loop1, 12) == 4 % delt del

        alpha = linedata(loop1, 13);
        zinfo = linedata(loop1:loop1 + 2, 3:5) + j * linedata(loop1:loop1 + 2, 6:8);
        yinfo = 1 / zinfo(1, 1);

        Y4 =- (1 / (3 * 1)) * [-2 * yinfo yinfo yinfo; yinfo -2 * yinfo yinfo; yinfo yinfo -2 * yinfo] + eye(3) * 0.00000;
        Y3 =- (1 / (3 * alpha * alpha)) * [-2 * yinfo yinfo yinfo; yinfo -2 * yinfo yinfo; yinfo yinfo -2 * yinfo] + eye(3) * 0.00000;
        Y1 =- (1 / (3 * alpha)) * [-2 * yinfo yinfo yinfo; yinfo -2 * yinfo yinfo; yinfo yinfo -2 * yinfo] + eye(3) * 0.00000;
        Y2 =- (1 / (3 * alpha)) * [-2 * yinfo yinfo yinfo; yinfo -2 * yinfo yinfo; yinfo yinfo -2 * yinfo] + eye(3) * 0.00000;

        Ynn(:, :, fromnode(loop1), tonode(loop1)) = Y3; % this sign needs be revised
        Ynn(:, :, tonode(loop1), fromnode(loop1)) = Y4;
        Ynm(:, :, fromnode(loop1), tonode(loop1)) = Y1;
        Ynm(:, :, tonode(loop1), fromnode(loop1)) = Y2; % this sign needs be revised

    end

    %% Del-Grounded Wye connection

    %% Grounded Wye-Del connection

    %% Ungrounded wye-Del connection
    if linedata(loop1, 12) == 3 % ungrounded wye-delta
        alpha = linedata(loop1, 13);
        zinfo = linedata(loop1:loop1 + 2, 3:5) + j * linedata(loop1:loop1 + 2, 6:8);
        % Check for singularity
        if rcond(zinfo) < 1e-12
            yinfo = 1 / zinfo(1, 1);
        else
            yinfo = zinfo \ eye(3);
        end

        % Matrix for ungrounded wye side (no zero sequence current path)
        Yw = (1 / (3 * alpha ^ 2)) * [2 -1 -1; -1 2 -1; -1 -1 2] * yinfo;

        % Matrix for delta side
        Yd = (1/3) * [2 -1 -1; -1 2 -1; -1 -1 2] * yinfo;

        % Mutual admittance matrices (connection transformation)
        Ywd = (1 / (3 * alpha)) * [2 -1 -1; -1 2 -1; -1 -1 2] * yinfo;
        Ydw = Ywd; % symmetrical for reciprocity

        % Assign to system admittance matrices
        Ynn(:, :, fromnode(loop1), tonode(loop1)) = Yw;
        Ynn(:, :, tonode(loop1), fromnode(loop1)) = Yd;
        Ynm(:, :, fromnode(loop1), tonode(loop1)) = Ywd;
        Ynm(:, :, tonode(loop1), fromnode(loop1)) = Ydw;
    end

end

%% Consoldiated Y-bus
nodes = unique([fromnode; tonode]);
nodes(isnan(nodes)) = [];

for loop1 = 1:max(nodes)
    loop1
    Ysum = 0;

    for loop2 = 1:max(nodes)
        loop2
        Ysum = Ysum + Ynn(:, :, loop1, loop2);

        if loop1 ~= loop2
            Ybusmat(loop1 * 3 - 2:loop1 * 3, loop2 * 3 - 2:loop2 * 3) = -Ynm(:, :, loop1, loop2); % this sign needs be revised
        end

    end

    Ybusmat(loop1 * 3 - 2:loop1 * 3, loop1 * 3 - 2:loop1 * 3) = Ysum;
end

%% Write intermediate excel files
save('Ybusmat.mat', 'Ybusmat');
%xlswrite('Y',[real(Ybusmat) imag(Ybusmat)]);
