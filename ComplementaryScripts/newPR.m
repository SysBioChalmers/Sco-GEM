function newPR(model)
% newCommit
%   Prepares a new pull request (PR) of the scoGEM model to the devel
%   branch, for submission to GitHub. In contrast to the master branch,
%   the devel and derived branches do not contain .mat and .xlsx files.
%   This function should be run from the ComplementaryScripts directory.
%
%   model   RAVEN model structure to be used in PR. If no model
%           structure is provided, the stored XML file will be imported
%           prior to exporting other files
%
%   Usage: newPR(model)
%
% Eduard Kerkhoven, 2018-09-15

%Check if in master:
currentBranch = git('rev-parse --abbrev-ref HEAD');
if strcmp(currentBranch,'master')
    error('ERROR: in master branch. For new releases, use newRelease.m')
end

if ~exist('model')
    %Load model:
    model = importModel('../ModelFiles/xml/scoGEM.xml');
end

%Save model
exportForGit(model,'scoGEM','../',{'txt', 'xml', 'yml'});
end