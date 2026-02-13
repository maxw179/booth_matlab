% runRepeatedAcq.m

outputFolder = 'test';
if ~exist(outputFolder,'dir'); mkdir(outputFolder); end

% --- Create controllers ---
slm = SLMController( ...
    'BoardNumber', 1, ...
    'BitDepth', 12, ...
    'InvPiOver128', inv_pi_over_128, ...      % <-- your existing constant
    'LUT1024', 'C:\Program Files\Meadowlark Optics\Blink OverDrive Plus\LUT Files\slm6748_at1300_1stOrder_031325_3.lut', ...
    'LUT512',  'C:\Program Files\Meadowlark Optics\Blink OverDrive Plus\LUT Files\512x512_linearVoltage.LUT');

slm.load();

% Optional: mask (pick one)
slm.setCircularMask('Center',[slm.Width/2, slm.Height/2], 'Radius', 400);
% slm.setMask(uint8(ones(slm.Width*slm.Height,1))); % full aperture

si = ScanImageController(hSI, 'PollPeriodSec',0.2, 'FallbackSize',[si_height si_width]);

% --- Repeated acquisitions ---
N = 20;
pauseAfterSLM = 0.2;

for k = 1:N
    % Example E-field (replace with your real pattern)
    E = ones(slm.Width*slm.Height, 1);

    slm.writeEfield(E);
    pause(pauseAfterSLM);

    img = si.grabFrame();
    saveTiff16(img, fullfile(outputFolder, sprintf('Img_%04d.tif', k)));
end

% --- Cleanup ---
delete(slm);   % calls SDK delete + unload libs
