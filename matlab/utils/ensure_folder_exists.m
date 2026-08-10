function ensure_folder_exists(path, folder_name)
%ENSURE_FOLDER_EXISTS Create a child folder when it is missing.

    full_path = fullfile(path, folder_name);
    
    if ~exist(full_path, 'dir')
        mkdir(full_path);
        fprintf("Folder '%s' created at '%s'\n", folder_name, path);
    else
        fprintf("Folder '%s' already exists at '%s'\n", folder_name, path);
    end
end
