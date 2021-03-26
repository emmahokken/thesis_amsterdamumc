file = '../../../../../data/projects/ahead/raw_data/Subcortex_0005_002_R02/2018_06_21/Su_2919/su_21062018_0953202_7_2_wip_wb_inv1_senseV4.raw';


% read in data 
r = MRecon(file);

% select specific data type 
% 1 is standard data, 3 is phase correction data, 5 is noise 
r.Parameter.Parameter2Read.typ = 1;
% r.Parameter.Parameter2Read.Update;
r.ReadData;

% perform corrections 
r.DcOffsetCorrection;
r.PDACorrection;
r.RandomPhaseCorrection;
r.MeasPhaseCorrection;
r.SortData;
r.GridData;

d = r.Data;

f = fft2(d);
size(d)
imshow(f(:,:,140,1,1,1,1))