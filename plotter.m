% Open file
fid = fopen('objective.txt');

T = 10000;
% Construct y
y = zeros(T,1);

for i = 1:T
    tline = fgets(fid);
    C = strsplit(tline);
    for j = 1:length(C)
        if ~isempty(str2num(C{j}))
            y(i) = str2double(C{j});
            break;
        end
    end
end

% Close file
fclose(fid);

% Plot data and save figure
figure;
plot((1:T)', y);
saveas(gcf, './figures/convergence_curve.png')
