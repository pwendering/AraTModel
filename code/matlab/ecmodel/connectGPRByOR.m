function newGPR = connectGPRByOR(gpr)
% Splits GPR rules by ') AND (', which can be problematic with GECKO
% Joins all possible combinations by OR.
% Requires Neural Network Toolbox Matlab toolbox.
% 
% Input:        char gpr:       original GPR rule
% Output:       char newGPR:    new GPR rule

% split the original rule by AND-connector
splitPattern = '\) \& [x(]';
splitGPR = regexp(gpr,splitPattern,'split');

% get gene numbers from all components
geneIdx = cellfun(@(x)str2double(regexp(x,'\d+','match')),splitGPR,'un',0);

% get all possible combinations between elements within these components
combMat = geneIdx{1};
for i=2:numel(geneIdx)
    combMat = combvec(combMat,geneIdx{i});
end
    
% create new rule
newGPR = cell(1,size(combMat,2));
for i=1:size(combMat,2)
    tmpGPRPart = strcat(strcat('x(',strtrim(cellstr(num2str(combMat(:,i))))),')');
    newGPR{i} = ['( ' strjoin(tmpGPRPart,' & ') ' )'];
end
newGPR = strjoin(newGPR,' | ');

end