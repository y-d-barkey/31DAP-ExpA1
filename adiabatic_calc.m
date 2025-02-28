clear

gas = 'helium';
gas = 'compressed air';

if strcmp(gas, 'helium')
    folders = {'1.0', '2.0', '2.4'};
else 
    folders = {'1.0', '2.0', '3.0'};
end

Complete_table = table;
p_table = table;

% Ambient pressure and temperature at time of testing
p_amb = 120555;
t_amb = 296.03;
p_amb_SD = 1000
t_amb_SD = 0.214;

% Conversion from V to Pa
v_to_pa = @ (v) p_amb + (v + 0.05063) / 4.49503E-6
vsd_to_pasd = @ (v, vsd) sqrt( ...
    (vsd/4.49503E-6)^2 ...
    +(0.03266/4.49503E-6)^2 ...
    +(2.14197E-7*(v + 0.05063)/(4.49503E-6)^2)^2 ...
    );

% Conversion from V to K
v_to_k = @ (v) 273.15 + (v-3.0588)/(0.0218881118881);
vsd_to_ksd = @ (v, vsd) sqrt(0.00583261377445+(vsd/0.0218881118881)^2);

% Parse through each folder (Initial pressure)
for j = 1:3
    % Change the directory if needed
    directory = strcat( ...
    'C:\31DAP\250218_ExpA\experiment A\pressure release measurements\', ...
    gas, ...
    '\', ...
    folders{j}, ...
    ' bar\' ...
    );
    
    files = dir(fullfile(directory, '*.txt'));
    Ad_const_table = table;
    
    
    % Calculate the Adiabatic constant (AB and BC) for each measurement
    for i = 1:length(files)
        filename = fullfile(directory, files(i).name);
        
        if ~files(i).isdir
            
            % Preparing the data
            adia_data = readtable(filename, "Delimiter", "\t");
            vps = adia_data.Voltage_A0TG;
            vts = adia_data.Voltage_A1TG;
            
            % Finding the time of pressure release
            pdropoff = find(diff(vps) < -0.015);
            prelease = pdropoff(1);
            
            % Finding the time of temperature release
            tdropoff = find(diff(vts) < -0.003);
            trelease = tdropoff(1);
            
            % Calculate mean and std.dev for voltage pressure at point A
            vpA_range = vps(1:prelease - 2);
            vpA = mean(vpA_range);
            vpA_SD = std(vpA_range);
            
            % Calculate mean and std.dev for voltage temp at point A
            vtA_range = vts(1:trelease - 2);
            vtA = mean(vtA_range);
            vtA_SD = std(vtA_range);
            
            % Select the bc trajectory out of the data
            Tt_bc = adia_data.Time_sec2(find(vts == min(vts)):end);
            vt_bc = vts(find(vts == min(vts)):end);
            Tp_bc = adia_data.Time_sec1(find(vts == min(vts)):end);
            vp_bc = vps(find(vts == min(vts)):end);
    
            % Fit the data to an exponential decay
            bc_fit_function = @(v, x) v(1) + (v(2)) * exp(-v(3) * x);
            vp0 = [vp_bc(end), -0.5, 0.2];
            vp_bc_fit_model = fitnlm(Tp_bc, vp_bc, bc_fit_function, vp0);
            % The fit for asymptotic temperature ended up unnecessary
            % vt0 = [vt_bc(end), -0.5, 0.2];
            % vt_bc_fit_model = fitnlm(Tp_bc, vt_bc, bc_fit_function, vt0);
            
            fit = bc_fit_function( ...
                vp_bc_fit_model.Coefficients.Estimate, Tp_bc ...
                );
            
            % Pressure at point A
            pA = v_to_pa(vpA);
            pA_SD = vsd_to_pasd(vpA, vpA_SD);
            
            % Temperature at point A
            tA = v_to_k(vtA);
            tA_SD = vsd_to_ksd(vpA, vtA_SD);
            
            % Temperature at point B
            tB = v_to_k(min(vts));
            tB_SD = vsd_to_ksd(min(vts), min(vts) * vtA_SD / vtA);
    
            % Pressure at point C
            pC = v_to_pa(vp_bc_fit_model.Coefficients.Estimate(1));
            pC_SD = vsd_to_pasd( ...
                vp_bc_fit_model.Coefficients.Estimate(1), ...
                vp_bc_fit_model.Coefficients.SE(1) ...
                );
            
            % Temperature at point C (PERHAPS NOT NEEDED?)
            % tC = v_to_k(vt_bc_fit_model.Coefficients.Estimate(1));
            % tC_SD = vsd_to_ksd( ...
                %vt_bc_fit_model.Coefficients.Estimate(1), ...
                %vt_bc_fit_model.Coefficients.SE(1) ...
                %);       
            
            % Adiabatic const. for trajectory AB (Adiabatic Expansion)
            g_ab = 1 / (1 - log(tB/tA) / log(p_amb/pA));
            g_ab_SD = g_ab^2/log(p_amb/pA) * sqrt( ...
                (tB_SD/tB)^2+ ...
                (tA_SD/tA)^2 + ...
                ((pA_SD/pA)^2+(p_amb_SD/p_amb)^2)* ...
                (log(tB/tA)/log(p_amb/pA))^2 ...
                );
            
            % Adiabatic const. for trajectory BC (Isochoric Heating)
            g_bc = log(pA/p_amb) / log(pA/pC);
            g_bc_SD = (1/log(pA/pC))*sqrt( ...
                (pA_SD*(1-g_bc)/pA)^2 + ...
                (pC_SD*g_bc/pC)^2 + ...
                (p_amb_SD/p_amb)^2 ...
                );
            
            cofaab = [ ...
                (tB_SD/tB)^2, ...
                (tA_SD/tA)^2, ...
                ((pA_SD/pA)*log(tB/tA)/log(p_amb/pA))^2, ...
                ((p_amb_SD/p_amb)*log(tB/tA)/log(p_amb/pA))^2 ...
            ];
            cofabc = [ ...
                (pA_SD*(1-g_bc)/pA)^2, ...
                (pC_SD*g_bc/pC)^2, ...
                (p_amb_SD/p_amb)^2 ...
            ];

            %disp('AB:');
            %disp(find(cofaab == max(cofaab)));
            %disp('BC:');
            %disp(find(cofabc == max(cofabc)));

            

            % Add values to the table
            newrow = table(g_ab, g_ab_SD, g_bc, g_bc_SD);
            Ad_const_table = [Ad_const_table; newrow];

            prow = table(str2double(folders{j}), pA, pA_SD, pC, pC_SD);
            p_table = [p_table; prow];
        end    
    end
    
    % Weighted average for AB
    g_ab_w = 1 ./ (Ad_const_table.g_ab_SD).^2;
    g_ab_wmean = sum(g_ab_w .* Ad_const_table.g_ab) / ...
        sum(g_ab_w);
    g_ab_wSD = sqrt( ...
        (6/5) * sum(g_ab_w .* (Ad_const_table.g_ab - g_ab_wmean).^2) / ...
        sum(g_ab_w) ...
        );
    
    % Weighted average for BC
    g_bc_w = 1 ./ Ad_const_table.g_bc_SD.^2;
    g_bc_wmean = sum(g_bc_w .* Ad_const_table.g_bc) / ...
        sum(g_bc_w);
    g_bc_wSD = sqrt( ...
        (6/5) * sum(g_bc_w .* (Ad_const_table.g_bc - g_bc_wmean).^2) / ...
        sum(g_bc_w) ...
        );
    
    % Add the result for the adiabatic constant to the table
    complete_row = table( ...
        g_ab_wmean, g_ab_wSD, ...
        g_bc_wmean, g_bc_wSD ...
        );

    Complete_table = [Complete_table; complete_row];
end

Complete_table = addvars(Complete_table, str2double(folders)', 'NewVariableNames', 'p_initial', 'Before', 'g_ab_wmean');
disp(Complete_table)
% disp(p_table)

% Prepare variables for graph
x = Complete_table.p_initial;
y1 = Complete_table.g_ab_wmean;
y1_err = Complete_table.g_ab_wSD;
y2 = Complete_table.g_bc_wmean;
y2_err = Complete_table.g_bc_wSD;

w1 = 1./(y1_err).^2;
w2 = 1./(y2_err).^2;

% Calculate weighted averages
g_ab_wm = sum(w1 .* y1) / sum(w1);
g_bc_wm = sum(w2 .* y2) / sum(w2);
g_wm = (sum(w1 .* y1) + sum(w2 .* y2) ) / (sum(w1) + sum(w2));



% Calculate weighted std.dev from weighted variances
var_ab = sum(w1 .* (y1 - g_wm).^2);
var_bc = sum(w2 .* (y2 - g_wm).^2);
g_wSD = sqrt((6/5)*(var_ab + var_bc) / (sum(w1) + sum(w2)));

% Calculate weighted std.dev individually
g_ab_wSD = sqrt((3/2)*sum(w1 .* (y1 - g_ab_wm).^2)/sum(w1));
g_bc_wSD = sqrt((3/2)*sum(w2 .* (y2 - g_bc_wm).^2)/sum(w2));

all_gamma = [y1; y2];
all_gamma_SD = [y1_err;y2_err];

disp(['Adiabatic constant AB: ', num2str(g_ab_wm),char(177),num2str(g_ab_wSD)])
disp(['Adiabatic constant BC: ', num2str(g_bc_wm),char(177),num2str(g_bc_wSD)])
disp(['Adiabatic constant: ', num2str(g_wm),char(177),num2str(g_wSD)])
disp(['Reduced Chi-Squared: ', num2str(sum(((all_gamma - g_wm)./all_gamma_SD).^2) / 5)])




hold on;
errorbar(x, y1, y1_err)
errorbar(x, y2, y2_err)
yline(g_wm+g_wSD)
yline(g_wm-g_wSD)
yline(g_wm)
hold off;