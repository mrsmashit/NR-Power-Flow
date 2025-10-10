clc;
clearvars;

Ybusmat = load('Ybusmat.mat');
Ybusmat = Ybusmat.Ybusmat;

Ymag = abs(Ybusmat);
Yang = angle(Ybusmat);
org_size_Y = length(Ymag);
all_zeros = find(all(Ymag == 0));
Ybusmat(:, all_zeros) = [];
Ybusmat(all_zeros, :) = [];

Ymag(:, all_zeros) = [];
Ymag(all_zeros, :) = [];
n = length(Ymag);
Yang(:, all_zeros) = [];
Yang(all_zeros, :) = [];

loaddata = xlsread('IEEE_4bus_data.xls', 'Load', 'A1:H2');

PloadA = [0; 0; 0; loaddata(1, 3)];
PloadB = [0; 0; 0; loaddata(1, 4)];
PloadC = [0; 0; 0; loaddata(1, 5)];
QloadA = [0; 0; 0; loaddata(1, 6)];
QloadB = [0; 0; 0; loaddata(1, 7)];
QloadC = [0; 0; 0; loaddata(1, 8)];

for i = 1:length(PloadA)
    Pload(3 * i - 2:3 * i) = [PloadA(i) PloadB(i) PloadC(i)];
    Qload(3 * i - 2:3 * i) = [QloadA(i) QloadB(i) QloadC(i)];
    Pgen(3 * i - 2:3 * i) = [0 0 0];
    Qgen(3 * i - 2:3 * i) = [0 0 0];
end

Pload = Pload' * 1e3; Qload = Qload' * 1e3; Pgen = Pgen' * 1e3; Qgen = Qgen' * 1e3;
Pload(all_zeros) = []; Qload(all_zeros) = [];
Pgen(all_zeros) = []; Qgen(all_zeros) = [];

% Finding which Jacobians are required
bust = ones(1, n) * 3;
bust(1:3) = 1;

busdata = xlsread('IEEE_4bus_data.xls', 'Bus', 'A1:C5');
Vnom_pri = busdata(2, 2);
Vnom_sec = busdata(4, 2);

Vfull = ones(1, org_size_Y);

Vfull(1:3:6) = Vnom_pri * 1e3 * Vfull(1:3:6) * exp(j * 0 * pi / 180);
Vfull(2:3:6) = Vnom_pri * 1e3 * Vfull(2:3:6) * exp(j * -120 * pi / 180);
Vfull(3:3:6) = Vnom_pri * 1e3 * Vfull(3:3:6) * exp(j * 120 * pi / 180);

% Include 30° phase shift for delta side
Vfull(7:3:12) = Vnom_sec * 1e3 * Vfull(7:3:12) * exp(j * 0 * pi / 180);
Vfull(8:3:12) = Vnom_sec * 1e3 * Vfull(8:3:12) * exp(j * (-120) * pi / 180);
Vfull(9:3:12) = Vnom_sec * 1e3 * Vfull(9:3:12) * exp(j * 120 * pi / 180);

V = abs(Vfull);
delta = angle(Vfull);
V(all_zeros) = []; delta(all_zeros) = [];

tic

for p = 1:n

    if bust(p) == 1
        M(p) = 0;
    else
        M(p) = 1;
    end

end

for p = n + 1:2 * n

    if bust(p - n) == 3
        M(p) = 1;
    else
        M(p) = 0;
    end

end

for p = 1:2 * n

    for q = 1:2 * n
        J(p, q) = M(p) * M(q);
    end

end

error = 100; % use while loop to run until error is less than specifeid
itcount = 0; %iteration count

while error > 0.001 % specified error

    if itcount == 500
        break;
    end

    % for certain number of iteration use for loop below, comment while
    % loop
    %for loop=1:10 % start of main loop

    itcount = itcount + 1;
    cnt = 1;
    % For some display
    errorv(itcount) = error;
    Vnode3(itcount) = V(3);

    % Finding the mismatch equations, 1 equation each for PV bus and 2
    % equations each for PQ bus
    for p = 1:n

        if ((bust(p) == 2) | (bust(p) == 3))
            sum1 = 0;

            for q = 1:n
                sum1 = sum1 + V(q) * Ymag(p, q) * cos(delta(p) - delta(q) - Yang(p, q));
            end

            Mismt(cnt) = Pgen(p) - Pload(p) - V(p) * sum1;
            cnt = cnt + 1;
        end

    end

    for p = 1:n

        if (bust(p) == 3)
            sum1 = 0;

            for q = 1:n
                sum1 = sum1 + V(q) * Ymag(p, q) * sin(delta(p) - delta(q) - Yang(p, q));
            end

            Mismt(cnt) = Qgen(p) - Qload(p) - V(p) * sum1;
            cnt = cnt + 1;
        end

    end

    % start of jacobaina calculation, stored in J1

    ctr1 = 1;

    for p = 1:n

        if (M(p) == 1)
            ctr2 = 1;

            for q = 1:n

                if (J(p, q) == 1)
                    sum2 = 0;

                    if (p == q)

                        for r = 1:n
                            sum2 = sum2 + V(r) * Ymag(p, r) * sin(delta(p) - delta(r) - Yang(p, r));
                        end

                        J1(ctr1, ctr2) = V(p) * sum2 + V(p) * V(p) * Ymag(p, p) * sin(Yang(p, p));
                    else
                        J1(ctr1, ctr2) = -V(p) * V(q) * Ymag(p, q) * sin(delta(p) - delta(q) - Yang(p, q));
                    end

                    ctr2 = ctr2 + 1;
                end

            end

            for q = n + 1:2 * n

                if (J(p, q) == 1)
                    sum2 = 0;

                    if (p == q - n)

                        for r = 1:n
                            sum2 = sum2 + V(r) * Ymag(p, r) * cos(delta(p) - delta(r) - Yang(p, r));
                        end

                        J1(ctr1, ctr2) = -sum2 - V(p) * Ymag(p, p) * cos(Yang(p, p));
                    else
                        J1(ctr1, ctr2) = -V(p) * Ymag(p, q - n) * cos(delta(p) - delta(q - n) - Yang(p, q - n));
                    end

                    ctr2 = ctr2 + 1;
                end

            end

            ctr1 = ctr1 + 1;

        end

    end

    for p = n + 1:2 * n

        if (M(p) == 1)
            ctr2 = 1;

            for q = 1:n

                if (J(p, q) == 1)
                    sum2 = 0;

                    if (p - n == q)

                        for r = 1:n
                            sum2 = sum2 + V(p - n) * V(r) * Ymag(p - n, r) * cos(delta(p - n) - delta(r) - Yang(p - n, r));
                        end

                        J1(ctr1, ctr2) = -sum2 + V(p - n) * V(p - n) * Ymag(p - n, p - n) * cos(Yang(p - n, p - n));
                    else
                        J1(ctr1, ctr2) = V(p - n) * V(q) * Ymag(p - n, q) * cos(delta(p - n) - delta(q) - Yang(p - n, q));
                    end

                    ctr2 = ctr2 + 1;
                end

            end

            for q = n + 1:2 * n

                if (J(p, q) == 1)
                    sum2 = 0;

                    if (p - n == q - n)

                        for r = 1:n
                            sum2 = sum2 + V(r) * Ymag(p - n, r) * sin(delta(p - n) - delta(r) - Yang(p - n, r));
                        end

                        J1(ctr1, ctr2) = -sum2 + V(p - n) * Ymag(p - n, p - n) * sin(Yang(p - n, p - n));
                    else
                        J1(ctr1, ctr2) = -V(p - n) * Ymag(p - n, q - n) * sin(delta(p - n) - delta(q - n) - Yang(p - n, q - n));
                    end

                    ctr2 = ctr2 + 1;
                end

            end

            ctr1 = ctr1 + 1;
        end

    end % End of jacobian calcualtion

    % matrix inversion to calculate  the increments
    Inc = -1 * inv(J1) * Mismt';
    % Combined vector of V and delta

    Var = [delta V];
    ctr3 = 0;

    for p = 1:2 * n

        if (M(p) == 1)
            ctr3 = ctr3 + 1;
            Var(p) = Var(p) + Inc(ctr3);
        end

    end

    delta = Var(1:n);

    V = Var(n + 1:2 * n);
    error = max(abs(Mismt));

end % end of main loop

tf = toc;

text(5, .93, strcat('Time Taken (NR)   =', num2str(tf)))

Vnonlinfull = ones(org_size_Y, 1);

for i = 1:org_size_Y
    Vnonlinfull(all_zeros) = 0;
    Vnonlinfull(setdiff(1:end, all_zeros)) = V;
end

deltafull = ones(org_size_Y, 1);

for i = 1:org_size_Y
    deltafull(all_zeros) = 0;
    deltafull(setdiff(1:end, all_zeros)) = delta;
end

VphsA = Vnonlinfull(1:3:end);
VphsB = Vnonlinfull(2:3:end);
VphsC = Vnonlinfull(3:3:end);

Vnonlinfull = Vnonlinfull .* exp(j * deltafull);

% Calculate transformer currents considering wye-delta configuration
% Primary side (wye) currents
Ipri =- (Ybusmat(4:6, 1:3)) * (Vnonlinfull(1:3) - Vnonlinfull(4:6));

% Secondary side (delta) currents - convert line currents to phase currents
Isec_line = (Ybusmat(7:9, 4:6)) * (Vnonlinfull(4:6) - Vnonlinfull(7:9));
% Convert line currents to phase currents for delta connection
Isec_phase = zeros(3, 1);
Isec_phase(1) = (Isec_line(1) - Isec_line(3)) / sqrt(3); % Phase a
Isec_phase(2) = (Isec_line(2) - Isec_line(1)) / sqrt(3); % Phase b
Isec_phase(3) = (Isec_line(3) - Isec_line(2)) / sqrt(3); % Phase c

% Calculate apparent power
Ssub_pri = Vnonlinfull(1:3) .* conj(Ipri);
Ssub_sec = Vnonlinfull(4:6) .* conj(Isec_line);

% Save bus voltages in per unit and angles in degrees
% Get the valid indices (non-zero buses)
valid_indices = setdiff(1:length(Vnonlinfull), all_zeros);

% Create bus numbers array (1,1,1,2,2,2,etc. for each phase)
bus_numbers = ceil(valid_indices / 3)';

% Create base voltage array
V_base = zeros(size(Vnonlinfull));
V_base(1:6) = Vnom_pri * 1e3; % Primary side voltage base
V_base(7:end) = Vnom_sec * 1e3; % Secondary side voltage base

% Calculate per unit values and angles
V_pu = (abs(Vnonlinfull(valid_indices)) ./ V_base(valid_indices));
V_angle_deg = angle(Vnonlinfull(valid_indices)) * 180 / pi;

if size(V_angle_deg, 1) == 1
    V_angle_deg = V_angle_deg';
end

% Add phase identifiers (1=A, 2=B, 3=C)
phase_ids = mod(valid_indices - 1, 3) + 1;

% Combine results into a matrix [bus_number phase_id voltage_pu angle_deg]
voltage_results = [bus_numbers phase_ids' V_pu V_angle_deg];

% Save to .mat file
save('voltage_results.mat', 'voltage_results');
%}
