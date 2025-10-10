clc;
clearvars;
%% Read Line Data

linedata = xlsread('IEEE_4bus_data.xls', 'Feeder', 'A1:D4');
fromnode = linedata(:, 1);
tonode = linedata(:, 2);
llen = linedata(:, 3) * 0.000189394; % in miles
lconf = linedata(:, 4); % configuration no.
config = xlsread('IEEE_4bus_data.xls', 'Config', 'A2:S3');

%% Feeder configuration Data Processing

confno = length(config(:, 1));

for loop = 1:confno

    zconf = [config(loop, 2:4) + j * config(loop, 8:10)
           0 config(loop, 5:6) + j * config(loop, 11:12)
           0 0 config(loop, 7) + j * config(loop, 13)];
    bconf = [j * config(loop, 14:16)
           0 j * config(loop, 17:18)
           0 0 j * config(loop, 19)] * 1e-6;
    zinfo(:, :, config(loop, 1)) = lowerfill(zconf);
    binfo(:, :, config(loop, 1)) = lowerfill(bconf);

end

%% ABCD of feeder sections
zlineinfo = cell(length(linedata), 14);
I3 = eye(3);
feederno = length(fromnode);

for loop1 = 1:length(fromnode) % loops over all connections

    for loop2 = 1:length(zinfo) % loops over all feeder configurations

        if lconf(loop1) == loop2

            zlineinfo(loop1 * 3 - 2, 1:5) = [num2cell([fromnode(loop1) tonode(loop1) real(zinfo(1, 1:3, loop2)) * llen(loop1)]')]';
            zlineinfo(loop1 * 3 - 1, 1:5) = [num2cell([fromnode(loop1) tonode(loop1) real(zinfo(2, 1:3, loop2)) * llen(loop1)]')]';
            zlineinfo(loop1 * 3, 1:5) = [num2cell([fromnode(loop1) tonode(loop1) real(zinfo(3, 1:3, loop2)) * llen(loop1)]')]';
            zlineinfo(loop1 * 3 - 2, 6:8) = [num2cell((imag(zinfo(1, 1:3, loop2)) * llen(loop1))')]';
            zlineinfo(loop1 * 3 - 1, 6:8) = [num2cell((imag(zinfo(2, 1:3, loop2)) * llen(loop1))')]';
            zlineinfo(loop1 * 3, 6:8) = [num2cell((imag(zinfo(3, 1:3, loop2)) * llen(loop1))')]';
            zlineinfo(loop1 * 3 - 2, 9:11) = [num2cell((imag(binfo(1, 1:3, loop2)) * llen(loop1))')]';
            zlineinfo(loop1 * 3 - 1, 9:11) = [num2cell((imag(binfo(2, 1:3, loop2)) * llen(loop1))')]';
            zlineinfo(loop1 * 3, 9:11) = [num2cell((imag(binfo(3, 1:3, loop2)) * llen(loop1))')]';
            zlineinfo{loop1 * 3 - 1, 1} = [];
            zlineinfo{loop1 * 3 - 1, 2} = [];
            zlineinfo{loop1 * 3, 1} = [];
            zlineinfo{loop1 * 3, 2} = [];

            zlineinfo((loop1) * 3 - 2, 12) = num2cell(1); %feeder

            break

        end

    end

end

%% Read Transformer Data
trdata = xlsread('IEEE_4bus_data.xls', 'Transformer', 'A2:P2');
l1 = size(trdata);

if (l1(1) >= 1) & (l1(2) >= 1) % if transformer exists

    fromnode = trdata(:, 1);
    tonode = trdata(:, 2);

    for loop1 = 1:length(fromnode) % loops over all transformer
        %% grounded Wye- Grounded Wye connection
        if trdata(loop1, 9) == 1 & trdata(loop1, 10) == 1
            ntr = [trdata(loop1, 7) / trdata(loop1, 8) trdata(loop1, 7) / trdata(loop1, 8) trdata(loop1, 7) / trdata(loop1, 8)]';

            Ztr = (trdata(loop1, 11) / 100) * (trdata(loop1, 8) * trdata(loop1, 8) / (trdata(loop1, 3) / 1000));
            Zti = (trdata(loop1, 12) / 100) * (trdata(loop1, 8) * trdata(loop1, 8) / (trdata(loop1, 3) / 1000));

            zlineinfo((feederno + loop1) * 3 - 2, 12) = num2cell(2); % transformer
            zlineinfo((feederno + loop1) * 3 - 2, 13) = num2cell(ntr(1)); % wye wye

        end

        %% Del-Del connection

        if trdata(loop1, 9) == 2 & trdata(loop1, 10) == 2
            ntr = trdata(loop1, 7) / trdata(loop1, 8);

            Ztr = (trdata(loop1, 11) / 100) * (trdata(loop1, 8) * trdata(loop1, 8) / (trdata(loop1, 3) / 1000)); % assumes idential impedence on each phase

            Zti = (trdata(loop1, 12) / 100) * (trdata(loop1, 8) * trdata(loop1, 8) / (trdata(loop1, 3) / 1000));

            zlineinfo((feederno + loop1) * 3 - 2, 12) = num2cell(4); % tr
            zlineinfo((feederno + loop1) * 3 - 2, 13) = num2cell(ntr(1)); % deldel

        end

        %% Del-Grounded Wye connection

        %% Grounded Wye-Del connection

        %% Ungrounded wye-Del connection
        if trdata(loop1, 9) == 3 && trdata(loop1, 10) == 3
            ntr = [trdata(loop1, 7) / trdata(loop1, 8) trdata(loop1, 7) / trdata(loop1, 8) trdata(loop1, 7) / trdata(loop1, 8)].' * exp(-1i * pi / 6);

            Zpct = (trdata(loop1, 11) + 1i * trdata(loop1, 12)) / 100;
            Zbase = (trdata(loop1, 8) ^ 2) / (trdata(loop1, 3) / 1000);
            Z = Zpct * Zbase;
            Ztr = real(Z);
            Zti = imag(Z);

            zlineinfo((feederno + loop1) * 3 - 2, 12) = num2cell(3); % transformer
            zlineinfo((feederno + loop1) * 3 - 2, 13) = num2cell(ntr(1)); % wye del
        end

        %% Assigning to zlineifo

        zlineinfo((feederno + loop1) * 3 - 2, 1:5) = [num2cell([fromnode(loop1) tonode(loop1) [Ztr 0 0]]')]';
        zlineinfo((feederno + loop1) * 3 - 1, 1:5) = [num2cell([fromnode(loop1) tonode(loop1) [0 Ztr 0]]')]';
        zlineinfo((feederno + loop1) * 3, 1:5) = [num2cell([fromnode(loop1) tonode(loop1) [0 0 Ztr]]')]';
        zlineinfo((feederno + loop1) * 3 - 2, 6:8) = [num2cell([Zti 0 0]')]';
        zlineinfo((feederno + loop1) * 3 - 1, 6:8) = [num2cell([0 Zti 0]')]';
        zlineinfo((feederno + loop1) * 3, 6:8) = [num2cell([0 0 Zti]')]';
        zlineinfo((feederno + loop1) * 3 - 2, 9:11) = [num2cell([0 0 0]')]';
        zlineinfo((feederno + loop1) * 3 - 1, 9:11) = [num2cell([0 0 0]')]';
        zlineinfo((feederno + loop1) * 3, 9:11) = [num2cell([0 0 0]')]';
        zlineinfo{(feederno + loop1) * 3 - 1, 1} = [];
        zlineinfo{(feederno + loop1) * 3 - 1, 2} = [];
        zlineinfo{(feederno + loop1) * 3, 1} = [];
        zlineinfo{(feederno + loop1) * 3, 2} = [];

    end

end

%% Write to excel files

xlswrite('zinfo_4bus.xls', zlineinfo);
