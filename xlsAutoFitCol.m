%% SOURCE COPYRIGHT

% Copyright (c) 2008, Tal Shir
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

function xlsAutoFitCol(filename,sheetname,varargin)

%set the column to auto fit

%xlsAutoFitCol(filename,sheetname,range)
% Example:
%xlsAutoFitCol('filename','Sheet1','A:F')

range = varargin{1};
    
[fpath,file,ext] = fileparts(char(filename));
if isempty(fpath)
    fpath = pwd;
end
Excel = actxserver('Excel.Application');
set(Excel,'Visible',0);
Workbook = invoke(Excel.Workbooks, 'open', [fpath filesep file ext]);

sheet = get(Excel.Worksheets, 'Item',sheetname);
invoke(sheet,'Activate');

ExAct = Excel.Activesheet;
   
ExActRange = get(ExAct,'Range',range);
ExActRange.Select;

invoke(Excel.Selection.Columns,'Autofit');

invoke(Workbook, 'Save');
invoke(Excel, 'Quit');
delete(Excel);

