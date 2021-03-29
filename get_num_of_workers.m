function get_num_of_workers

p = gcp('nocreate');
if isempty(p)
    n_workers = 0;
else
    n_workers = p.NumWorkers;
end
disp(['Number of workers currently used: ', num2str(n_workers)]);

warning('off', 'MATLAB:DELETE:Permission');
myCluster = parcluster('local');
warning('on', 'MATLAB:DELETE:Permission');
disp(['Total number of workers available: ', num2str(myCluster.NumWorkers)]);

