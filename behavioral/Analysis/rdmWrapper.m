function RDM = rdmWrapper(mat, cond, subID, color)
% mat: trial x trial dissimilarity  matrix
% cond: string for condition
% subID: string for subject ID
% color: 1x3 rgb color
RDM.RDM = mat;
RDM.name = [cond ' | ' subID ' | Session: 1'];
if nargin < 4
    RDM.color = [0 0 1];
else
    RDM.color = color;

end