baphy_set_path
dbopen;

sql = ['SELECT gSingleCell.*,gPenetration.animal'...
  ' FROM gSingleCell,gPenetration'...
  ' WHERE gSingleCell.penid=gPenetration.id'...
  ' AND area in ("A1") AND animal in ("Tango") ORDER BY cellid'];
%pfc="PFC", "FC", "RPFC", "LPFC"
%a1: "A1"  %, "AC", "LAC", "RAC", "LA1", "RA1"
%peg="PEG", "PPF", "PSF", "border"
%VPr= "proPPF"
celldata=mysql(sql);

cell_list = {celldata.cellid};
runclass = '"CLK"'; % CCH: can be reference only and reference/target. Leave only ref/tar.

for i = 1 : length(cell_list)
  fprintf('\n#########################                        cell %i of %i\n', i, length(cell_list));
   A = findrecfromcell(cell_list{1,i}, runclass, 0);
   if ~isempty(A)
     cell_list{2,i} = A;
   end
end


fprintf('\n\nFound %i cellids\n', length(cell_list));



%% delete columns with no recording info
row1 = cell_list(1,:);
row2 = cell_list(2,:);

row1 = row1(~cellfun(@isempty, row2)); % delete empty elements
row2 = row2(~cellfun(@isempty, row2));

cell_list = [row1;row2];
clear row1;
clear row2;

fprintf('Only %i cellids did task %s\n', length(cell_list), runclass);

%% get rid of CCHs with no target (complex chord as reference only, instead of torc vs tone)

if strcmp(runclass, 'CCH')
row1 = cell_list(1,:);
row2 = cell_list(2,:);


  for i = 1 : length(row2)
    discardCCH = [];
    for j = 1 : size(row2{i},2)
      if isempty(row2{i}{j}) || ...
          ~strcmp(row2{i}{j}.exptparams.TrialObject.ReferenceHandle.descriptor, 'Torc') || ...
          isempty(row2{i}{j}.exptparams.TrialObject.TargetHandle)
        discardCCH = [discardCCH, j];
      end
      
    end
    for del = fliplr(discardCCH)
      row2{i}(del) = [];
    end
  end
  cell_list = [row1;row2];
end
clear row1;
clear row2;

row1 = cell_list(1,:);
row2 = cell_list(2,:);

row1 = row1(~cellfun(@isempty, row2)); % delete empty elements
row2 = row2(~cellfun(@isempty, row2));

cell_list = [row1;row2];

clear row1;
clear row2;

% %% now leave only those with any active behavior
% row1 = cell_list(1,:);
% row2 = cell_list(2,:);
% 
% idxanyactive = []; % indices with any active behavior
% for i = 1 : length(row2)
%   found = 0; % used to avoid repeating cellids with more than one active instance
%   for j = 1: length(row2{i})
%     if ~isempty(row2{i}{j}) && strcmp(row2{i}{j}.Behavior, 'a') && found == 0
%       idxanyactive = [idxanyactive i];
%       found = 1;
%     end
%   end
% end
% 
% cell_list = [row1(idxanyactive);row2(idxanyactive)];
% fprintf('%i cellids were tested on active behavior\n\n', length(cell_list));
% 
% clear row1; clear row2;


%% now save passive-active sequence
row1 = cell_list(1,:);
row2 = cell_list(2,:);

row3 = cell(1,length(row2));

for i = 1 : length(row2)
  row3{i} = char(zeros(1,length(row2{i})));
  for j = 1 : length(row2{i})
    row3{i}(j) = row2{i}{j}.Behavior;
  end
end

cell_list = [cell_list; row3];
