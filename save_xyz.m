function save_xyz(fname, my_matrix)
    fid = fopen(fname,'wt');
    for i = 1:size(my_matrix,1)
        fprintf(fid,'%g ',my_matrix(i,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
end