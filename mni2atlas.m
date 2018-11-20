function [atlas]=mni2atlas(roi,atlas_selector,thr_or_res)
%MNI2ATLAS: takes in input a ROI or a vector of coordinates (both in the MNI
%space) and returns labels from different fsl atlases.
%
%In VECTOR modality labels are returned in probability values (same results
%of fsl atlas tool).
%In ROI modality the probability value reported for a label represents the
%frequency of that label in the roi for a given threshold of fsl atlas
%probability map.
%__________________________________________________________________________
%HOW TO USE
%   MNI2ATLAS(VECTOR/ROI) the first input can be a MNI vector or a ROI in
%   the MNI space. Depending on the input the script switch between two
%   different work modalities. With no other input the script will seek labels
%   among all accesible fsl atalses.
%
%   MNI2ATLAS(VECTOR/ROI,ATLAS_SELECTOR) allow to choose between the 
%   following atales:
%       1) Juelich Histological Atlas
%       2) Harvard-Oxford Cortical Structural Atlas
%       3) Harvard-Oxford Subcortical Structural Atlas
%       4) JHU ICBM-DTI-81 White Matter labels
%       5) JHU White Matter tractography Atlas
%       6) Oxford Thalamic Connectivity Atlas
%       7) Cerebellar Atlas in MNI152 after FLIRT
%       8) Cerebellar Atlas in MNI152 after FNIRT
%   ATLAS_SELECTOR must be a raw vector (i.e. [1,3,6]). Default value is
%   [1:1:8]. You can also leave it as an empty vector (i.e. (VECTOR/ROI,[])).
%
%   [ATLAS]=MNI2ATLAS(VECTOR/ROI,...) the script returns the structure          
%   ATLAS whit the following fields: .name (of the atlas), .labels (a cell
%   vector). No stout will be print.
%
%   MNI2ATLAS(VECTOR) prints on screen labels found for the mni VECTOR
%   position.
%
%   MNI2ATLAS(ROI) prints on screen labels found for the input ROI. ROI can 
%   be a preloaded (with load_nii) volume or the path of a nii volume. 
%
%ADVANCED OPTIONS
%   MNI2ATLAS(ROI,ATLAS_SELECTOR,THR) THR allows to choose among three 
%   different threshold levels: 0, 25, 50 (e.i. 0%, 25%, 50%). Default
%   value is 25. Option available only under ROI modality.
%
%   MNI2ATLAS(VECTOR,ATLAS_SELECTOR,RESOLUTION) RESOLUTION allows to
%   choose between 1mm or 2mm atlases. 1mm atlases performe better regions
%   identification but requires more loading time. Default value is 2mm.
%   Option available only under VECTOR modality.
%__________________________________________________________________________
%SYSTEM REQUIREMENTS
%   1) run only on linux system
%   2) NifTI and ANALYZE tool (version > 2012-10-12)
%   3) fsl atlases
%__________________________________________________________________________
%
%   Version 1.0
%   Last modified 11/10/2013
%   danielemascali@gmail.com

%end of help


global roi2atlas_atlas_roi_struct    roi2atlas_roi_check_input        %for ROI
global roi2atlas_atlas_vector_struct roi2atlas_vector_check_input  %for vector

%--------------------------------------------------------------------------
%                      system/tool check
%--------------------------------------------------------------------------

if exist('load_nii') == 0 
    fprintf('>This function requires NifTI/ANALYZE tool.\n You can download the last version from:\n http://research.baycrest.org/~jimmy/NIfTI/');
    return
end

if ~isunix
    fprintf('>This function was written for linux system only.\n If you really want to use a different system you can upgrade the script replacing\n the linux system function "locate" (in the nasted functions: load_atlas_vector\n and load_atlas_roi)  with an equivalent one for your system.\n');
    return    
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%                      argin check
%--------------------------------------------------------------------------

if nargin == 0
    help mni2atlas
    return
end

if nargin < 2 || isempty(atlas_selector)
    
    atlas_selector = 1:1:8;

end

if nargin < 3 || isempty(thr_or_res)
    
    thr = 25;
    resolution = '2mm';  %you can also put '1mm' but the result will be almost the same but much more time demanding
    
elseif ischar(thr_or_res)
    
    if strcmp(thr_or_res,'1mm') || strcmp(thr_or_res,'2mm')
    
        resolution = thr_or_res;
        
    else
        disp('>No valid resolution. Please enter the string 1mm or 2mm.')
        return
    end
    
else
    
    if thr_or_res == 0 || thr_or_res == 25 || thr_or_res == 50
    
        thr = thr_or_res;
       
    else
        disp('>No valid threshold. Please enter one of the following values: 0, 25, 50.')
        return
    end    
    
end

if ischar(roi) || (~ischar(roi) && numel(roi) ~= 3)
    
    if nargin >= 3
        
        if ischar(thr_or_res)
            disp('>In roi modality the third input must be a scalar (0, 20 or 50) not a string.')
            return
        end
        
    end
    
elseif ~ischar(roi) && numel(roi) == 3
    
    if nargin >= 3
        
        if ~ischar(thr_or_res)
           disp('>In vector modality the third input must be a string (1mm or 2mm).')
           return
        end
        
    end
    
end


if ischar(roi)
    
    if isempty(dir(roi))
        
        disp('>The input roi is not in the current or specified path.')
        return
        
    end
    
    roi = load_nii(roi);
    roi = uint8(roi.img);
    
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


print_info.modality = [];       %print out struct
print_info.atlas_selector = atlas_selector;
print_info.roi.thr = [];
print_info.vector.resolution = [];



%--------------------------------------------------------------------------
%                       BODY 
%--------------------------------------------------------------------------
if  numel(roi) == 3 % è stato inserito un vettore mni non una roi
        
    print_info.modality = 'vec';
    print_info.vector.resolution = resolution;
    
    if isempty(roi2atlas_vector_check_input) || (~isempty(roi2atlas_vector_check_input) && (~strcmp(roi2atlas_vector_check_input.resolution,resolution)))
        
       
        error = load_atlas_vector(resolution);  
        
        if error == 1
            return
        end

        if ~isfield(roi2atlas_atlas_vector_struct,'volume')
            
            disp('>No fsl atlases found in your system.');
            clearvars -global roi2atlas_vector_check_input roi2atlas_vector_check_input
            return
        
        end
        
    end

    
    choosen_atlas_struct = select_atlas_vector(atlas_selector);

    [xyz_cord] = mni2xyz(roi,str2double(resolution(1)));
    
    output_struct = vector_process(xyz_cord,choosen_atlas_struct);
    
elseif ndims(roi) == 3 % è stata inserita una roi
    
    
    %-------check roi resolution----------------   
    if sum(size(roi) == [91,109,91]) == 3
        
        resolution = '2mm';
        
    elseif sum(size(roi) == [182,218,182]) == 3
        
        resolution = '1mm';
        
    else

        disp('>The size of you ROI is not compatible with the geometry of fsl.\nCompatible geometries are [91,109,91] and [182,218,182].');
        return

    end
    %-------------------------------------------
    
    print_info.modality = 'roi';
    print_info.roi.thr = thr;
    
    if isempty(roi2atlas_roi_check_input) || (~isempty(roi2atlas_roi_check_input) && (roi2atlas_roi_check_input.thr ~= thr || ~strcmp(roi2atlas_roi_check_input.resolution,resolution)))
               
        error = load_atlas_roi(thr,resolution);
        
        if error == 1
            return
        end
        
        if ~isfield(roi2atlas_atlas_roi_struct,'volume')
            
            disp('>No fsl atlases found in your system');
            clearvars -global roi2atlas_roi_check_input roi2atlas_roi_check_input
            return
        
        end
        
    end
        
    choosen_atlas_struct = select_atlas_roi(atlas_selector);
    
    output_struct = roi_process(roi,choosen_atlas_struct);
    
else
   
    disp('>The first input is neither a roi neither a MNI position.')
    return
    
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
%                       PRINT/OUTPUT
%--------------------------------------------------------------------------
if nargout == 1
    
    atlas = prepare_output(output_struct);      %no stout
  
else
    
    print_labels(output_struct,print_info)  %no variable storage, yes stout
    
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

return 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           ROI FUNCTIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [error] = load_atlas_roi(thr,resolution)

global roi2atlas_atlas_roi_struct roi2atlas_roi_check_input

error =0;

roi2atlas_atlas_roi_struct(1).private_name = ['Juelich-maxprob-thr',num2str(thr),'-',resolution,'.nii.gz'];
%roi2atlas_atlas_roi_struct(2).private_name = ['HarvardOxford-cortl-maxprob-thr',num2str(thr),'-',resolution,'.nii.gz'];
roi2atlas_atlas_roi_struct(2).private_name = ['HarvardOxford-cort-maxprob-thr',num2str(thr),'-',resolution,'.nii.gz'];
roi2atlas_atlas_roi_struct(3).private_name = ['HarvardOxford-sub-maxprob-thr',num2str(thr),'-',resolution,'.nii.gz '];
roi2atlas_atlas_roi_struct(4).private_name = ['JHU-ICBM-labels-',resolution,'.nii.gz'];
roi2atlas_atlas_roi_struct(5).private_name = ['JHU-ICBM-tracts-maxprob-thr',num2str(thr),'-',resolution,'.nii.gz'];
roi2atlas_atlas_roi_struct(6).private_name = ['Thalamus-maxprob-thr',num2str(thr),'-',resolution,'.nii.gz'];
roi2atlas_atlas_roi_struct(7).private_name = ['Cerebellum-MNIflirt-maxprob-thr',num2str(thr),'-',resolution,'.nii.gz'];
roi2atlas_atlas_roi_struct(8).private_name = ['Cerebellum-MNIfnirt-maxprob-thr',num2str(thr),'-',resolution,'.nii.gz'];
%roi2atlas_atlas_roi_struct(9).private_name = ['Talairach-labels-',resolution,'.nii.gz'];


roi2atlas_atlas_roi_struct(1).name = 'Juelich Histological Atlas' ;
roi2atlas_atlas_roi_struct(2).name = 'Harvard-Oxford Cortical Structural Atlas';
roi2atlas_atlas_roi_struct(3).name = 'Harvard-Oxford Subcortical Structural Atlas';
roi2atlas_atlas_roi_struct(4).name = 'JHU ICBM-DTI-81 White Matter labels';
roi2atlas_atlas_roi_struct(5).name = 'JHU White Matter tractography Atlas';
roi2atlas_atlas_roi_struct(6).name = 'Oxford Thalamic Connectivity Atlas';
roi2atlas_atlas_roi_struct(7).name = 'Cerebellar Atlas in MNI152 after FLIRT';
roi2atlas_atlas_roi_struct(8).name = 'Cerebellar Atlas in MNI152 after FNIRT';
%roi2atlas_atlas_roi_struct(9).name = 'Talairach Daemon Atlas';

roi2atlas_atlas_roi_struct(1).xml_name = 'Juelich.xml';
%roi2atlas_atlas_roi_struct(2).xml_name = 'HarvardOxford-Cortical-Lateralized.xml';
roi2atlas_atlas_roi_struct(2).xml_name = 'HarvardOxford-Cortical.xml';
roi2atlas_atlas_roi_struct(3).xml_name = 'HarvardOxford-Subcortical.xml';
roi2atlas_atlas_roi_struct(4).xml_name = 'JHU-labels.xml';
roi2atlas_atlas_roi_struct(5).xml_name = 'JHU-tracts.xml';
roi2atlas_atlas_roi_struct(6).xml_name = 'Thalamus.xml';
roi2atlas_atlas_roi_struct(7).xml_name = 'Cerebellum_MNIflirt.xml';
roi2atlas_atlas_roi_struct(8).xml_name = 'Cerebellum_MNIfnirt.xml';
%roi2atlas_atlas_roi_struct(9).xml_name = 'Talairach.xml';

tot_step = length(roi2atlas_atlas_roi_struct);
h = waitbar(0,'Loading atlases, be patient...');

for i=1:length(roi2atlas_atlas_roi_struct)
    
    waitbar(i/tot_step);
    
    [null,result] = system(['locate ',roi2atlas_atlas_roi_struct(i).private_name]); % le vecchie versioni di matlab non riconoscono la tilde, ho messo null come varibile di non interesse per non creare problemi di compatibilit�.
    
    if isempty(result)
        
        warning(['Can not find atlas: ',roi2atlas_atlas_roi_struct(i).private_name]);
        
    else
    
        path = result(1:(strfind(result,'.nii.gz') + 6));
    
        try
            tmp = load_nii(path);
        catch
            fprintf('>Seems you have an old version of NifTI/ANALYZE tool that is not able to open .gz file.\n You can download the last version from:\n http://research.baycrest.org/~jimmy/NIfTI/\n');
            close(h)
            error = 1;
            return
        end
        
    
        roi2atlas_atlas_roi_struct(i).volume = tmp.img;
    
       
    
        [null,result] = system(['locate ',roi2atlas_atlas_roi_struct(i).xml_name]);
        
        if isempty(result)
            
            warning(['Can not find atlas labels: ',roi2atlas_atlas_roi_struct(i).xml_name]);
            
        else

            path = result(1:(strfind(result,'.xml') + 3));

            xDoc =xmlread(path);

            roi2atlas_atlas_roi_struct(i).xml_loaded = xmlwrite(xDoc);
            
        end
    
    
    end
end

close(h);

roi2atlas_roi_check_input.thr = thr;
roi2atlas_roi_check_input.resolution = resolution;

return
end

function [a] = select_atlas_roi(atlas_selector)

global roi2atlas_atlas_roi_struct

a=roi2atlas_atlas_roi_struct;

for i = length(a):-1:1
    
    if sum(i==atlas_selector) == 0
        
        a(i) = [];
                
    end
    
end

return
end

function [a] = roi_process(roi,a)

for i=1:length(a)
    
    temp = a(i).volume(roi == 1);
    
    label_index = unique(temp);
    label_freq = nan(length(label_index),1);
    
    total_voxel = length(temp);
    
    for j=1:length(label_index)
        
        if label_index(j) == 0
            continue
        end
        
        count = sum(sum(temp == label_index(j)));
        
        label_freq(j) = count/total_voxel;
        
    end
    
    
    label_freq(label_index == 0) = [];
    label_index(label_index == 0) = [];
         
    a(i).label = [];   
    
    [sorted_label_freq,sorted_label_freq_index] = sort(label_freq,'descend');
    
    for j=1:length(sorted_label_freq)
        
        if isempty(label_index)
            
            a(i).label{j} = [];
            break
            
        end
               
        start_index = strfind(a(i).xml_loaded,['index="',num2str(label_index(sorted_label_freq_index(j)) - 1),'"']);
        
        tmp = a(i).xml_loaded(start_index:end);
        
        tmp_inf = strfind(tmp,'>');
        tmp_sup = strfind(tmp,'<');
        
        a(i).label{j} = [num2str(sorted_label_freq(j)*100,'%2.0f'),'% ',tmp(tmp_inf(1) +1 :tmp_sup(1) - 1)];
        
           
    end
    
        
end

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           VECTOR FUNCTIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [error] = load_atlas_vector(resolution)

global roi2atlas_atlas_vector_struct roi2atlas_vector_check_input

error = 0;

roi2atlas_atlas_vector_struct(1).private_name = ['Juelich-prob-',resolution,'.nii.gz'];
roi2atlas_atlas_vector_struct(2).private_name = ['HarvardOxford-cortl-prob-',resolution,'.nii.gz'];
roi2atlas_atlas_vector_struct(3).private_name = ['HarvardOxford-sub-prob-',resolution,'.nii.gz '];
roi2atlas_atlas_vector_struct(4).private_name = ['JHU-ICBM-labels-',resolution,'.nii.gz'];
roi2atlas_atlas_vector_struct(5).private_name = ['JHU-ICBM-tracts-prob-',resolution,'.nii.gz'];
roi2atlas_atlas_vector_struct(6).private_name = ['Thalamus-prob-',resolution,'.nii.gz'];
roi2atlas_atlas_vector_struct(7).private_name = ['Cerebellum-MNIflirt-prob-',resolution,'.nii.gz'];
roi2atlas_atlas_vector_struct(8).private_name = ['Cerebellum-MNIfnirt-prob-',resolution,'.nii.gz'];
%roi2atlas_atlas_vector_struct(9).private_name = ['Talairach-labels-',resolution,'.nii.gz'];


roi2atlas_atlas_vector_struct(1).name = 'Juelich Histological Atlas' ;
roi2atlas_atlas_vector_struct(2).name = 'Harvard-Oxford Cortical Structural Atlas';
roi2atlas_atlas_vector_struct(3).name = 'Harvard-Oxford Subcortical Structural Atlas';
roi2atlas_atlas_vector_struct(4).name = 'JHU ICBM-DTI-81 White Matter labels';
roi2atlas_atlas_vector_struct(5).name = 'JHU White Matter tractography Atlas';
roi2atlas_atlas_vector_struct(6).name = 'Oxford Thalamic Connectivity Atlas';
roi2atlas_atlas_vector_struct(7).name = 'Cerebellar Atlas in MNI152 after FLIRT';
roi2atlas_atlas_vector_struct(8).name = 'Cerebellar Atlas in MNI152 after FNIRT';
%roi2atlas_atlas_vector_struct(9).name = 'Talairach Daemon Atlas';

roi2atlas_atlas_vector_struct(1).xml_name = 'Juelich.xml';
roi2atlas_atlas_vector_struct(2).xml_name = 'HarvardOxford-Cortical-Lateralized.xml';
roi2atlas_atlas_vector_struct(3).xml_name = 'HarvardOxford-Subcortical.xml';
roi2atlas_atlas_vector_struct(4).xml_name = 'JHU-labels.xml';
roi2atlas_atlas_vector_struct(5).xml_name = 'JHU-tracts.xml';
roi2atlas_atlas_vector_struct(6).xml_name = 'Thalamus.xml';
roi2atlas_atlas_vector_struct(7).xml_name = 'Cerebellum_MNIflirt.xml';
roi2atlas_atlas_vector_struct(8).xml_name = 'Cerebellum_MNIfnirt.xml';
%roi2atlas_atlas_vector_struct(9).xml_name = 'Talairach.xml';

tot_step = length(roi2atlas_atlas_vector_struct);
h = waitbar(0,'Loading atlases, be patient...');

for i=1:length(roi2atlas_atlas_vector_struct)
    
    waitbar(i/tot_step);
    
    [null,result] = system(['locate ',roi2atlas_atlas_vector_struct(i).private_name]);
    
    if isempty(result)
        
        warning(['Can not find atlas: ',roi2atlas_atlas_vector_struct(i).private_name]);
        
    else
    
        path = result(1:(strfind(result,'.nii.gz') + 6));
    
        try
            tmp = load_nii(path);
        catch
            fprintf('>Seems you have an old version of NifTI/ANALYZE tool that is not able to open .gz file\n You can download the last version from:\n http://research.baycrest.org/~jimmy/NIfTI/\n');
            close(h)
            error = 1;
            return
        end
        
        roi2atlas_atlas_vector_struct(i).volume = tmp.img;
    
       
    
        [null,result] = system(['locate ',roi2atlas_atlas_vector_struct(i).xml_name]);
        
        if isempty(result)
            
            warning(['Can not find atlas labels: ',roi2atlas_atlas_vector_struct(i).xml_name]);
            
        else

            path = result(1:(strfind(result,'.xml') + 3));

            xDoc =xmlread(path);

            roi2atlas_atlas_vector_struct(i).xml_loaded = xmlwrite(xDoc);
            
        end
    
    
    end

end

close(h);

roi2atlas_vector_check_input.resolution = resolution;

return
end

function [a] = select_atlas_vector(atlas_selector)

global roi2atlas_atlas_vector_struct

a=roi2atlas_atlas_vector_struct;

for i = length(a):-1:1
    
    if sum(i==atlas_selector) == 0
        
        a(i) = [];
                
    end
    
end

return
end

function [a] = vector_process(xyz,a)

for i=1:length(a)
    
     
    temp = a(i).volume(xyz(1),xyz(2),xyz(3),:);
    
    temp =squeeze(temp);
    
    a(i).label = [];   
    
    if sum(temp) == 0
        
        continue 
        
    else
        
        count = 0;
        
        [temp_sort,temp_index] = sort(temp,'descend');
        
        for j=1:length(temp)
            
            if temp_sort(j) == 0
                break
            end
          
            count = count +1;

            start_index = strfind(a(i).xml_loaded,['index="',num2str(temp_index(j) - 1),'"']);

            tmp = a(i).xml_loaded(start_index:end);

            tmp_inf = strfind(tmp,'>');
            tmp_sup = strfind(tmp,'<');

            a(i).label{count} = [num2str(temp_sort(j),'%2.0f'),'% ',tmp(tmp_inf(1) +1 :tmp_sup(1) - 1)];
                
            
        end
    end
end
                
return
end

function [xyz] = mni2xyz(mni,vs)

if vs == 2

    origin = [45 63 36]; % [X Y Z]
    
    mni(1) = vs*round(mni(1)/vs);
    mni(2) = vs*fix(mni(2)/vs);
    mni(3) = vs*fix(mni(3)/vs);

elseif vs == 1
    
    origin = [91 126 72]; % [X Y Z]

end

xyz(1)=origin(1) + mni(1)/vs +1;      %original was origin(1) - mni(1)/vs
xyz(2)=mni(2)/vs + origin(2) +1;
xyz(3)=mni(3)/vs + origin(3) +1;
            
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           OUTPUT FUNCTIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output] = prepare_output(a)

output = [];

for i=1:length(a)
    
    output(i).name = a(i).name;
    output(i).label = a(i).label;

end

return
end

function print_labels(a,b)

fprintf('_________________________________');
fprintf('\nMNI2ATLAS');
fprintf('\n>modality:\t\t');
if strcmp(b.modality,'vec')
    fprintf('vector');
    fprintf('\n>resolution:\t\t%s',b.vector.resolution)
else
    fprintf('roi');
    fprintf('\n>thrshold:\t\t%d%%',b.roi.thr)
end
fprintf('\n>atlas selected:\t%d\n',length(b.atlas_selector))
fprintf('_________________________________\n');
total_label_found = 0;

for i=1:length(a)
    
    if ~isempty(a(i).label)
        
        total_label_found = total_label_found + 1;
        
    end
    
end

if total_label_found == 0
    
    fprintf('\n>No labels found \n');
else
    fprintf('\n>Found labels in %d atlases:\n',total_label_found);
end

for i=1:length(a)
    
    if ~isempty(a(i).label) 
          
        fprintf('\n%s',a(i).name)
        
        for j=1:length(a(i).label)
        
            fprintf('\n\t%s',a(i).label{j})
            
        end
        
        fprintf('\n');
        
    end
end
        
fprintf('\n');
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
