function startup()
    run_startup();
    dirs = {'autodiff', 'model-io', 'multiscale', 'modules', 'solvers', ...
        'visualization','co2lab'};
    base_path = fileparts(mfilename('fullpath')); %#ok
    clear fn
    for i = 1:numel(dirs)
        mrstPath('addroot', fullfile(base_path, dirs{i}));
    end
    run_local();
end

function d = rootdir()
   d = fileparts(mfilename('fullpath'));
end

function run_startup()
   run(fullfile(rootdir(), 'core', 'startup.m'));
end

function run_local()
   do_run_local(fullfile(rootdir(), 'startup_user.m'));
end

function do_run_local(local)
   if exist(local, 'file') == 2
      run(local);
   end
end
