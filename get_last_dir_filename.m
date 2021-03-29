function last_dir_filename = get_last_dir_filename(fname)

if ~exist('fname', 'var')
    fname = dbstack;
    if ismember(fname(2).name, {'get_last_dir', 'get_last_file', 'get_file_name', 'get_dir_name'})
        fname = fname(3).name;
    else
        fname = fname(2).name;
    end
end

if ispc
    homedir = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
else
    homedir = getenv('HOME');
end

[~, me] = system('whoami');
me = me(find(me == '\', 1, 'last') + 1 : end-1);
last_dir_filename = fullfile(homedir, [fname, '_', me, '.mat']);
