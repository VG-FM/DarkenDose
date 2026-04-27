% developed by Vinicius Gomes of Brazil
% DarkenDose - Radiochromic film dose reconstruction in cGy

clear; clc; close all;

%% ===================== SETTINGS =====================
channelName   = 'red';     % 'red', 'green', 'blue', 'gray'
showCropOnly  = true;
doseCaxis     = [];        % e.g. [0 450], or [] for automatic
isoLineWidth  = 2.0;
labelSpacing  = 300;

exportCSV = true;
exportMAT = true;
exportPNG = true;
exportTXT = true;

doseUnit = 'cGy';

% Safe output name
t = clock;
outputBaseName = sprintf('DarkenDose_%04d%02d%02d_%02d%02d%02.0f', ...
    t(1), t(2), t(3), t(4), t(5), floor(t(6)));

% Calibration curve in cGy:
% OD(D) = (0.00128 + 0.00147*D) / (1 + 0.00165*D)
%
% Inverted curve:
% D(OD) = (0.00128 - OD) / (0.00165*OD - 0.00147)
%
% IMPORTANT:
% D is in cGy throughout the code.

%% ===================== INPUT FILES =====================
[doseFileName, dosePath] = uigetfile({'*.tif;*.tiff','TIFF Files (*.tif, *.tiff)'}, ...
    'Select irradiated film TIFF');
if isequal(doseFileName,0), error('Selection cancelled.'); end
doseFile = fullfile(dosePath, doseFileName);

[refFileName, refPath] = uigetfile({'*.tif;*.tiff','TIFF Files (*.tif, *.tiff)'}, ...
    'Select 0 cGy reference TIFF');
if isequal(refFileName,0), error('Selection cancelled.'); end
refFile = fullfile(refPath, refFileName);

fprintf('Dose file: %s\n', doseFile);
fprintf('Reference file: %s\n', refFile);

%% ===================== READ IMAGES =====================
doseInfo = imfinfo(doseFile);
refInfo  = imfinfo(refFile);

fprintf('\n=== DOSE TIFF INFO ===\n');
fprintf('ColorType: %s | BitDepth: %d | Size: %dx%d\n', ...
    doseInfo.ColorType, doseInfo.BitDepth, doseInfo.Width, doseInfo.Height);

fprintf('\n=== REFERENCE TIFF INFO ===\n');
fprintf('ColorType: %s | BitDepth: %d | Size: %dx%d\n', ...
    refInfo.ColorType, refInfo.BitDepth, refInfo.Width, refInfo.Height);

ImDose = double(imread(doseFile));
ImRef  = double(imread(refFile));

%% ===================== EXTRACT CHANNEL =====================
doseChannel = extractFilmChannel(ImDose, channelName);
refChannel  = extractFilmChannel(ImRef,  channelName);

doseChannel(doseChannel <= 0) = eps;
refChannel(refChannel <= 0)   = eps;

%% ===================== REFERENCE ROI =====================
figure('Name','Reference ROI');
imagesc(refChannel);
axis image;
colormap gray;
colorbar;
title('Click the top-left corner of the 0 cGy ROI');
drawnow;
[xa, ya] = ginput(1);

title('Click the bottom-right corner of the 0 cGy ROI');
drawnow;
[xb, yb] = ginput(1);

x1 = max(1, round(min(xa, xb)));
x2 = min(size(refChannel,2), round(max(xa, xb)));
y1 = max(1, round(min(ya, yb)));
y2 = min(size(refChannel,1), round(max(ya, yb)));

if x2 <= x1 || y2 <= y1
    error('Invalid ROI in reference image.');
end

refROI = refChannel(y1:y2, x1:x2);
I0 = mean(refROI(:), 'omitnan');

fprintf('\nI0 = %.6f\n', I0);

hold on;
plot([x1 x2 x2 x1 x1], [y1 y1 y2 y2 y1], 'y-', 'LineWidth', 2);
hold off;

%% ===================== FILM MASK =====================
figure('Name','Film Mask');
imagesc(doseChannel);
axis image;
colormap gray;
colorbar;
title('Click points around the film and press ENTER');
drawnow;
[xp, yp] = ginput;

if numel(xp) < 3
    error('At least 3 points are required to define a polygon.');
end

xp(end+1) = xp(1);
yp(end+1) = yp(1);

hold on;
plot(xp, yp, 'y-', 'LineWidth', 2);
hold off;

[XX, YY] = meshgrid(1:size(doseChannel,2), 1:size(doseChannel,1));
mask = inpolygon(XX, YY, xp, yp);

if ~any(mask(:))
    error('Empty mask. Please redraw the film contour.');
end

%% ===================== OPTICAL DENSITY =====================
OD = log10(I0 ./ doseChannel);
OD(~mask) = NaN;

%% ===================== DOSE CONVERSION =====================
denominator = 0.00165 .* OD - 0.00147;

Dose = (0.00128 - OD) ./ denominator;

Dose(abs(denominator) < 1e-8) = NaN;
Dose(Dose < 0) = NaN;
Dose(~mask) = NaN;

% Safety filter near singularity
OD_limit = 0.00147 / 0.00165;
Dose(OD >= OD_limit) = NaN;

%% ===================== DOSE STATISTICS =====================
doseValues = Dose(mask);
doseValues = doseValues(~isnan(doseValues));

if isempty(doseValues)
    error('No valid dose values were found inside the mask.');
end

meanDose   = mean(doseValues);
stdDose    = std(doseValues);
minDose    = min(doseValues);
maxDose    = max(doseValues);
medianDose = median(doseValues);

fprintf('\n=== DOSE STATISTICS ===\n');
fprintf('Mean dose    = %.4f %s\n', meanDose, doseUnit);
fprintf('Std dose     = %.4f %s\n', stdDose, doseUnit);
fprintf('Min dose     = %.4f %s\n', minDose, doseUnit);
fprintf('Max dose     = %.4f %s\n', maxDose, doseUnit);
fprintf('Median dose  = %.4f %s\n', medianDose, doseUnit);
fprintf('Valid pixels = %d\n', numel(doseValues));

%% ===================== CROP TO FILM REGION =====================
if showCropOnly
    [yy, xx] = find(mask);
    ymin = min(yy); ymax = max(yy);
    xmin = min(xx); xmax = max(xx);

    DosePlot = Dose(ymin:ymax, xmin:xmax);
else
    DosePlot = Dose;
end

%% ===================== ISODOSE LEVELS =====================
percentOffsets = -50:5:50;
isoLevels = meanDose * (1 + percentOffsets/100);

validIdx = isoLevels > 0;
isoLevels = isoLevels(validIdx);
percentOffsets = percentOffsets(validIdx);

fprintf('\n=== ISODOSE LEVELS BASED ON MEAN DOSE ===\n');
for i = 1:numel(isoLevels)
    fprintf('Isodose %2d: %8.4f %s (%+d%% from mean)\n', ...
        i, isoLevels(i), doseUnit, percentOffsets(i));
end

%% ===================== FIGURE 1: DOSE MAP =====================
fig1 = figure('Name','Dose Map');
imagesc(DosePlot);
axis image;
set(gca,'YDir','normal');
colormap(parula);
cb1 = colorbar;
ylabel(cb1, sprintf('Dose (%s)', doseUnit));
title(sprintf('Dose Map (%s)', doseUnit));

if ~isempty(doseCaxis)
    caxis(doseCaxis);
end

%% ===================== FIGURE 2: BLACK ISODOSES =====================
fig2 = figure('Name','Dose Map with Black Isodoses');
imagesc(DosePlot);
axis image;
set(gca,'YDir','normal');
colormap(parula);
cb2 = colorbar;
ylabel(cb2, sprintf('Dose (%s)', doseUnit));
title(sprintf('Dose Map with Isodose Contours - Mean = %.4f %s', meanDose, doseUnit));
hold on;

[C, hContour] = contour(DosePlot, isoLevels, ...
    'LineWidth', isoLineWidth, ...
    'LineColor', 'k');

clabel(C, hContour, ...
    'FontSize', 10, ...
    'Color', 'w', ...
    'FontWeight', 'bold', ...
    'LabelSpacing', labelSpacing);

hold off;

%% ===================== FIGURE 3: COLORED ISODOSES =====================
fig3 = figure('Name','Dose Map with Colored Isodoses');
imagesc(DosePlot);
axis image;
set(gca,'YDir','normal');
colormap(parula);
cb3 = colorbar;
ylabel(cb3, sprintf('Dose (%s)', doseUnit));
title(sprintf('Colored Isodose Contours - Mean = %.4f %s', meanDose, doseUnit));
hold on;

colors = turbo(numel(isoLevels));
legendHandles = gobjects(numel(isoLevels),1);
legendLabels  = strings(numel(isoLevels),1);

for i = 1:numel(isoLevels)
    contour(DosePlot, [isoLevels(i) isoLevels(i)], ...
        'LineWidth', isoLineWidth, ...
        'LineColor', colors(i,:));

    legendHandles(i) = plot(NaN, NaN, '-', ...
        'LineWidth', isoLineWidth, ...
        'Color', colors(i,:));

    legendLabels(i) = sprintf('%.4f %s (%+d%%)', ...
        isoLevels(i), doseUnit, percentOffsets(i));
end

legend(legendHandles, legendLabels, 'Location', 'eastoutside');
hold off;

%% ===================== EXPORT RESULTS =====================
if exportCSV
    writematrix(Dose, [outputBaseName '_DoseMap_cGy.csv']);
    writematrix(OD,   [outputBaseName '_ODMap.csv']);
    writematrix([percentOffsets(:), isoLevels(:)], ...
        [outputBaseName '_IsodoseLevels_cGy.csv']);
end

if exportMAT
    save([outputBaseName '_Results.mat'], ...
        'Dose', 'OD', 'mask', 'isoLevels', ...
        'I0', 'meanDose', 'stdDose', ...
        'minDose', 'maxDose', 'medianDose', ...
        'percentOffsets', 'doseUnit');
end

if exportPNG
    exportgraphics(fig1, [outputBaseName '_DoseMap_cGy.png'], 'Resolution', 300);
    exportgraphics(fig2, [outputBaseName '_DoseMap_BlackIsodoses_cGy.png'], 'Resolution', 300);
    exportgraphics(fig3, [outputBaseName '_DoseMap_ColoredIsodoses_cGy.png'], 'Resolution', 300);
end

if exportTXT
    fid = fopen([outputBaseName '_Report.txt'], 'w');

    fprintf(fid, 'DarkenDose Film Dose Report\n');
    fprintf(fid, 'Generated: %s\n\n', datestr(now));
    fprintf(fid, 'Channel: %s\n', channelName);
    fprintf(fid, 'Dose unit: %s\n\n', doseUnit);

    fprintf(fid, 'Calibration curve:\n');
    fprintf(fid, 'OD(D) = (0.00128 + 0.00147*D) / (1 + 0.00165*D)\n');
    fprintf(fid, 'D(OD) = (0.00128 - OD) / (0.00165*OD - 0.00147)\n');
    fprintf(fid, 'D is expressed in cGy.\n\n');

    fprintf(fid, 'I0: %.6f\n\n', I0);
    fprintf(fid, 'Mean dose: %.6f %s\n', meanDose, doseUnit);
    fprintf(fid, 'Std dose: %.6f %s\n', stdDose, doseUnit);
    fprintf(fid, 'Min dose: %.6f %s\n', minDose, doseUnit);
    fprintf(fid, 'Max dose: %.6f %s\n', maxDose, doseUnit);
    fprintf(fid, 'Median dose: %.6f %s\n', medianDose, doseUnit);
    fprintf(fid, 'Valid pixels: %d\n\n', numel(doseValues));

    fprintf(fid, 'Isodose levels based on mean dose:\n');
    fprintf(fid, 'Offset(%%),Dose(%s)\n', doseUnit);

    for i = 1:numel(isoLevels)
        fprintf(fid, '%+d,%.6f\n', percentOffsets(i), isoLevels(i));
    end

    fclose(fid);
end

disp('Processing completed.');

%% ===================== LOCAL FUNCTION =====================
function channel = extractFilmChannel(Im, channelName)

    if ndims(Im) == 2
        channel = Im;
        return;
    end

    switch lower(channelName)
        case 'red'
            channel = Im(:,:,1);
        case 'green'
            channel = Im(:,:,2);
        case 'blue'
            channel = Im(:,:,3);
        case 'gray'
            channel = mean(Im(:,:,1:3), 3);
        otherwise
            error('Invalid channel. Use red, green, blue, or gray.');
    end

end