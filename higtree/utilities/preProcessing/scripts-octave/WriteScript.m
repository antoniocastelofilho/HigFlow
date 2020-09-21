function WriteScript(DIM, ND, NB, NP, delta, reynolds, filedomain, fileboundary, fileprobe, filescript)

    filenameaux = sprintf('%s.sh',filescript);
    file        = fopen(filenameaux,'w');

    fprintf(file,'#!/bin/bash -f\n\n');
    fprintf(file,'amrs=$1\n');
    fprintf(file,'steps=$2\n');
    fprintf(file,'res=$3\n');
    fprintf(file,'mul=$4\n');
    fprintf(file,'vtks=$5\n\n');

    if DIM == 2
        fprintf(file,'\t./read-amr-solve-ns-example-2d \\\n');
    else
        if DIM == 3
            fprintf(file,'\t./read-amr-solve-ns-example-3d \\\n');
        end
    end

    fprintf(file,'\t\t%d \\\n',ND);
    for i = 1:ND
        fprintf(file,'\t\t\tamrs/%s-%d.amr \\\n',filedomain,i-1);
    end

    fprintf(file,'\t\t%d \\\n',NB);
    for i = 1:NB
        fprintf(file,'\t\t\tamrs/%s-%d.amr \\\n',fileboundary,i-1);
        if DIM == 2
            fprintf(file,'\t\t\t1 0 0 0 0 0 \\\n');
        else
            if DIM == 3
                fprintf(file,'\t\t\t1 0 0 0 0 0 0 0 \\\n');
            end
        end

        %fprintf(file,'\t\t\t\t1 0 0 0 0 0 \\\n');
    end

    fprintf(file,'\t\t%f %f $steps $res \\\n',delta,reynolds);

    fprintf(file,'\t\t%d \\\n',NP);
    for i = 1:NP
        fprintf(file,'\t\t\tamrs/%s-%d.amr \\\n',fileprobe,i-1);
    end

    fprintf(file,'\t\t$mul \\\n');
    fprintf(file,'\t\t$vtks \\\n');
    fprintf(file,'\t\t-ksp_type bcgs -pc_type hypre\n');

    fclose(file);

    return
end
