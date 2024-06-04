function subfolder_all = search_sub_folder(path_search_all)

dirinfo = dir(path_search_all);
dirinfo(~[dirinfo.isdir]) = [];
delete_idx = false(1, length(dirinfo));
for i = 1:length(dirinfo)
   if strcmp(dirinfo(i).name, '.')
       delete_idx(i) = true;
   elseif strcmp(dirinfo(i).name, '..')
       delete_idx(i) = true;
   end
end
dirinfo(delete_idx) = [];

subfolder_all = {};
if ~isempty(dirinfo)
    for i = 1:length(dirinfo)
        subfolder_temp = fullfile(dirinfo(i).folder, dirinfo(i).name);
        subfolder_all = [subfolder_all; subfolder_temp];
    end
end

end