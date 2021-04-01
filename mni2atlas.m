function [atlas]=mni2atlas(roi,atlas_selector,thr_or_res)
%MNI2ATLAS: from MNI VECTOR/ROI to FSL anatomical labels
%
%In VECTOR modality labels are returned in probability values (same results
%of FSL atlas tool).
%In ROI modality the probability value reported for a label represents the
%frequency of that label in the ROI for a given threshold of the FSL atlas
%probability map (0%, 25% or 50%; default = 25%).
%__________________________________________________________________________
%HOW TO USE
%   MNI2ATLAS(VECTOR/ROI) the first input can be a MNI vector or an ROI in
%   the MNI space. Depending on the input the script switches between two
%   different work modalities. With no other input the script will seek labels
%   among the available fsl atalses.
%
%   MNI2ATLAS(VECTOR/ROI,ATLAS_SELECTOR) allows to choose among the 
%   following atlases:
%       1) Juelich Histological Atlas
%       2) Harvard-Oxford Cortical Structural Atlas
%       3) Harvard-Oxford Subcortical Structural Atlas
%       4) JHU ICBM-DTI-81 White Matter labels
%       5) JHU White Matter tractography Atlas
%       6) Oxford Thalamic Connectivity Atlas
%       7) Cerebellar Atlas in MNI152 after FLIRT
%       8) Cerebellar Atlas in MNI152 after FNIRT
%       9) MNI Structural Atlas
%   ATLAS_SELECTOR must be a row vector (i.e. [1,3,6]). Default value is
%   [1:1:9]. You can also leave it as an empty vector (i.e., (VECTOR/ROI,[])).
%
%   [ATLAS]=MNI2ATLAS(VECTOR/ROI,...) the script returns the structure          
%   ATLAS with the following fields: .name (of the atlas), .labels (a cell
%   vector). No stdout will be print.
%
%   MNI2ATLAS(VECTOR) prints on screen labels found for the MNI VECTOR
%   position.
%
%   MNI2ATLAS(ROI) prints on screen labels found for the input ROI. ROI can 
%   be a preloaded (with load_nii) volume or the path to a nifti volume. 
%
%ADVANCED OPTIONS
%   MNI2ATLAS(ROI,ATLAS_SELECTOR,THR) THR allows to choose among 3 threshold
%   levels: 0, 25, 50 (i.e., 0%, 25%, 50%). Default
%   value is 25. Option available only under ROI modality.
%
%   MNI2ATLAS(VECTOR,ATLAS_SELECTOR,RESOLUTION) RESOLUTION allows to
%   choose between '1mm' or '2mm' atlases. 1mm atlases performs better region
%   identification but requires more loading time. Default value is '1mm'.
%   Option available only under VECTOR modality.
%__________________________________________________________________________
%SYSTEM REQUIREMENTS
%  <a href="matlab:
%web('https://it.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image')">NifTI and ANALYZE tool</a> (version > 2012-10-12)
%__________________________________________________________________________
%ACKNOWLEDGEMENTS
%  This function uses some of the available <a href="matlab:
%web('https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Atlases')">FSL atlases</a>.
%__________________________________________________________________________ 
%Daniele Mascali @ Enrico Fermi Center, MARBILab, Rome
%danielemascali@gmail.com

% First version: 2013
% Major update: 2018

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
    atlas_selector = 1:1:9;
end

if nargin < 3 || isempty(thr_or_res)
    %default values:
    thr = 25;
    resolution = '1mm';  
elseif ischar(thr_or_res)
    if strcmp(thr_or_res,'1mm') || strcmp(thr_or_res,'2mm')
        resolution = thr_or_res;
    else
        disp('>No valid resolution. Please enter the string ''1mm'' or ''2mm''.')
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
if  numel(roi) == 3 % the input is a vector (not an ROI) i.e., VECT mode
    print_info.modality = 'vec';
    print_info.vector.resolution = resolution;
    
    if isempty(roi2atlas_vector_check_input) || (~isempty(roi2atlas_vector_check_input) && (~strcmp(roi2atlas_vector_check_input.resolution,resolution)))
        load_atlas_vector(resolution);  
    end
    
    choosen_atlas_struct = select_atlas_vector(atlas_selector);
    [xyz_cord] = mni2xyz(roi,str2double(resolution(1)));
    output_struct = vector_process(xyz_cord,choosen_atlas_struct);
    
elseif ndims(roi) == 3 % the input is a ROI i.e., ROI mode
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
        load_atlas_roi(thr,resolution);
    end
        
    choosen_atlas_struct = select_atlas_roi(atlas_selector);
    output_struct = roi_process(roi,choosen_atlas_struct);
    
else
    disp('>The first input is neither a roi nor a MNI position.')
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
function load_atlas_roi(thr,resolution)

global roi2atlas_atlas_roi_struct roi2atlas_roi_check_input

[p,~,~] = fileparts(which('mni2atlas'));
base_path = [p,'/atlases/'];
roi2atlas_atlas_roi_struct(1).path = [base_path,'Juelich/Juelich-maxprob-thr',num2str(thr),'-',resolution,'.nii.gz'];
roi2atlas_atlas_roi_struct(2).path = [base_path,'HarvardOxford/HarvardOxford-cort-maxprob-thr',num2str(thr),'-',resolution,'.nii.gz'];
roi2atlas_atlas_roi_struct(3).path = [base_path,'HarvardOxford/HarvardOxford-sub-maxprob-thr',num2str(thr),'-',resolution,'.nii.gz'];
roi2atlas_atlas_roi_struct(4).path = [base_path,'JHU/JHU-ICBM-labels-',resolution,'.nii.gz'];
roi2atlas_atlas_roi_struct(5).path = [base_path,'JHU/JHU-ICBM-tracts-maxprob-thr',num2str(thr),'-',resolution,'.nii.gz'];
roi2atlas_atlas_roi_struct(6).path = [base_path,'Thalamus/Thalamus-maxprob-thr',num2str(thr),'-',resolution,'.nii.gz'];
roi2atlas_atlas_roi_struct(7).path = [base_path,'Cerebellum/Cerebellum-MNIflirt-maxprob-thr',num2str(thr),'-',resolution,'.nii.gz'];
roi2atlas_atlas_roi_struct(8).path = [base_path,'Cerebellum/Cerebellum-MNIfnirt-maxprob-thr',num2str(thr),'-',resolution,'.nii.gz'];
roi2atlas_atlas_roi_struct(9).path = [base_path,'MNI/MNI-maxprob-thr',num2str(thr),'-',resolution,'.nii.gz'];

roi2atlas_atlas_roi_struct(1).name = 'Juelich Histological Atlas' ;
roi2atlas_atlas_roi_struct(2).name = 'Harvard-Oxford Cortical Structural Atlas';
roi2atlas_atlas_roi_struct(3).name = 'Harvard-Oxford Subcortical Structural Atlas';
roi2atlas_atlas_roi_struct(4).name = 'JHU ICBM-DTI-81 White Matter labels';
roi2atlas_atlas_roi_struct(5).name = 'JHU White Matter tractography Atlas';
roi2atlas_atlas_roi_struct(6).name = 'Oxford Thalamic Connectivity Atlas';
roi2atlas_atlas_roi_struct(7).name = 'Cerebellar Atlas in MNI152 after FLIRT';
roi2atlas_atlas_roi_struct(8).name = 'Cerebellar Atlas in MNI152 after FNIRT';
roi2atlas_atlas_roi_struct(9).name = 'MNI Structural Atlas';

roi2atlas_atlas_roi_struct(1).xml_path = [base_path,'Juelich.xml'];
roi2atlas_atlas_roi_struct(2).xml_path = [base_path,'HarvardOxford-Cortical.xml'];
roi2atlas_atlas_roi_struct(3).xml_path = [base_path,'HarvardOxford-Subcortical.xml'];
roi2atlas_atlas_roi_struct(4).xml_path = [base_path,'JHU-labels.xml'];
roi2atlas_atlas_roi_struct(5).xml_path = [base_path,'JHU-tracts.xml'];
roi2atlas_atlas_roi_struct(6).xml_path = [base_path,'Thalamus.xml'];
roi2atlas_atlas_roi_struct(7).xml_path = [base_path,'Cerebellum_MNIflirt.xml'];
roi2atlas_atlas_roi_struct(8).xml_path = [base_path,'Cerebellum_MNIfnirt.xml'];
roi2atlas_atlas_roi_struct(9).xml_path = [base_path,'MNI.xml'];

tot_step = length(roi2atlas_atlas_roi_struct);
h = waitbar(0,'','name','Loading atlases, be patient...');

for l=1:length(roi2atlas_atlas_roi_struct)
    
    waitbar(l/tot_step,h,roi2atlas_atlas_roi_struct(l).name);
    
    tmp = load_nii(roi2atlas_atlas_roi_struct(l).path);
    
    roi2atlas_atlas_roi_struct(l).volume = tmp.img;
    
    xDoc =xmlread(roi2atlas_atlas_roi_struct(l).xml_path);
    
    roi2atlas_atlas_roi_struct(l).xml_loaded = xmlwrite(xDoc);
    
end

close(h);

roi2atlas_roi_check_input.thr = thr;
roi2atlas_roi_check_input.resolution = resolution;

return
end

function [a] = select_atlas_roi(atlas_selector)

global roi2atlas_atlas_roi_struct

a=roi2atlas_atlas_roi_struct;

for l = length(a):-1:1
    
    if sum(l==atlas_selector) == 0
        
        a(l) = [];
                
    end
    
end

return
end

function [a] = roi_process(roi,a)

for l=1:length(a)
    
    temp = a(l).volume(roi == 1);
    
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
         
    a(l).label = [];   
    
    [sorted_label_freq,sorted_label_freq_index] = sort(label_freq,'descend');
    
    for j=1:length(sorted_label_freq)
        
        if isempty(label_index)
            
            a(l).label{j} = [];
            break
            
        end
               
        start_index = strfind(a(l).xml_loaded,['index="',num2str(label_index(sorted_label_freq_index(j)) - 1),'"']);
        
        tmp = a(l).xml_loaded(start_index:end);
        
        tmp_inf = strfind(tmp,'>');
        tmp_sup = strfind(tmp,'<');
        
        a(l).label{j} = [num2str(sorted_label_freq(j)*100,'%2.0f'),'% ',tmp(tmp_inf(1) +1 :tmp_sup(1) - 1)];
        
           
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

function load_atlas_vector(resolution)

global roi2atlas_atlas_vector_struct roi2atlas_vector_check_input

[p,~,~] = fileparts(which('mni2atlas'));
base_path = [p,'/atlases/'];
roi2atlas_atlas_vector_struct(1).path = [base_path,'Juelich/Juelich-prob-',resolution,'.nii.gz'];
roi2atlas_atlas_vector_struct(2).path = [base_path,'HarvardOxford/HarvardOxford-cort-prob-',resolution,'.nii.gz'];
roi2atlas_atlas_vector_struct(3).path = [base_path,'HarvardOxford/HarvardOxford-sub-prob-',resolution,'.nii.gz'];
roi2atlas_atlas_vector_struct(4).path = [base_path,'JHU/JHU-ICBM-labels-',resolution,'.nii.gz'];
roi2atlas_atlas_vector_struct(5).path = [base_path,'JHU/JHU-ICBM-tracts-prob-',resolution,'.nii.gz'];
roi2atlas_atlas_vector_struct(6).path = [base_path,'Thalamus/Thalamus-prob-',resolution,'.nii.gz'];
roi2atlas_atlas_vector_struct(7).path = [base_path,'Cerebellum/Cerebellum-MNIflirt-prob-',resolution,'.nii.gz'];
roi2atlas_atlas_vector_struct(8).path = [base_path,'Cerebellum/Cerebellum-MNIfnirt-prob-',resolution,'.nii.gz'];
roi2atlas_atlas_vector_struct(9).path = [base_path,'MNI/MNI-prob-',resolution,'.nii.gz'];

roi2atlas_atlas_vector_struct(1).name = 'Juelich Histological Atlas' ;
roi2atlas_atlas_vector_struct(2).name = 'Harvard-Oxford Cortical Structural Atlas';
roi2atlas_atlas_vector_struct(3).name = 'Harvard-Oxford Subcortical Structural Atlas';
roi2atlas_atlas_vector_struct(4).name = 'JHU ICBM-DTI-81 White Matter labels';
roi2atlas_atlas_vector_struct(5).name = 'JHU White Matter tractography Atlas';
roi2atlas_atlas_vector_struct(6).name = 'Oxford Thalamic Connectivity Atlas';
roi2atlas_atlas_vector_struct(7).name = 'Cerebellar Atlas in MNI152 after FLIRT';
roi2atlas_atlas_vector_struct(8).name = 'Cerebellar Atlas in MNI152 after FNIRT';
roi2atlas_atlas_vector_struct(9).name = 'MNI Structural Atlas';

roi2atlas_atlas_vector_struct(1).xml_path = [base_path,'Juelich.xml'];
roi2atlas_atlas_vector_struct(2).xml_path = [base_path,'HarvardOxford-Cortical.xml'];
roi2atlas_atlas_vector_struct(3).xml_path = [base_path,'HarvardOxford-Subcortical.xml'];
roi2atlas_atlas_vector_struct(4).xml_path = [base_path,'JHU-labels.xml'];
roi2atlas_atlas_vector_struct(5).xml_path = [base_path,'JHU-tracts.xml'];
roi2atlas_atlas_vector_struct(6).xml_path = [base_path,'Thalamus.xml'];
roi2atlas_atlas_vector_struct(7).xml_path = [base_path,'Cerebellum_MNIflirt.xml'];
roi2atlas_atlas_vector_struct(8).xml_path = [base_path,'Cerebellum_MNIfnirt.xml'];
roi2atlas_atlas_vector_struct(9).xml_path = [base_path,'MNI.xml'];

tot_step = length(roi2atlas_atlas_vector_struct);
h = waitbar(0,'','name','Loading atlases, be patient...');

for l=1:length(roi2atlas_atlas_vector_struct)
    
    waitbar(l/tot_step,h,roi2atlas_atlas_vector_struct(l).name);
    
    tmp = load_nii(roi2atlas_atlas_vector_struct(l).path);
    
    roi2atlas_atlas_vector_struct(l).volume = tmp.img;
    
    xDoc =xmlread(roi2atlas_atlas_vector_struct(l).xml_path);

    roi2atlas_atlas_vector_struct(l).xml_loaded = xmlwrite(xDoc);

end

close(h);

roi2atlas_vector_check_input.resolution = resolution;

return
end

function [a] = select_atlas_vector(atlas_selector)

global roi2atlas_atlas_vector_struct

a=roi2atlas_atlas_vector_struct;

for l = length(a):-1:1
    if sum(l==atlas_selector) == 0  
        a(l) = [];      
    end
end

return
end

function [a] = vector_process(xyz,a)

for l=1:length(a)
    
    temp = a(l).volume(xyz(1),xyz(2),xyz(3),:);
    
    temp =squeeze(temp);
    
    a(l).label = [];   
    
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

            start_index = strfind(a(l).xml_loaded,['index="',num2str(temp_index(j) - 1),'"']);

            tmp = a(l).xml_loaded(start_index:end);

            tmp_inf = strfind(tmp,'>');
            tmp_sup = strfind(tmp,'<');

            a(l).label{count} = [num2str(temp_sort(j),'%2.1f'),'% ',tmp(tmp_inf(1) +1 :tmp_sup(1) - 1)];
                
        end
    end
end
                
return
end

function [xyz] = mni2xyz(mni,vs)

if vs == 2
    origin = [45 63 36]; % [X Y Z]
elseif vs == 1
    origin = [91 126 72]; % [X Y Z]
end


xyz(1)=origin(1) + round(mni(1)/vs) +1;      %was origin(1) - mni(1)/vs
xyz(2)=origin(2) + round(mni(2)/vs) +1;
xyz(3)=origin(3) + round(mni(3)/vs) +1;
            
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

for l=1:length(a)
    output(l).name = a(l).name;
    output(l).label = a(l).label;
end

return
end

function print_labels(a,b)

fprintf('___________________________________________');
fprintf('\n\t\tMNI2ATLAS');
fprintf('\n> Modality:\t\t');
if strcmp(b.modality,'vec')
    fprintf('vector');
    fprintf('\n> Resolution:\t\t%s',b.vector.resolution)
else
    fprintf('roi');
    fprintf('\n> Threshold:\t\t%d%%',b.roi.thr)
end
fprintf('\n> Mumber of atlases:\t%d\n',length(b.atlas_selector))
fprintf('___________________________________________\n');
total_label_found = 0;

for l=1:length(a)
    if ~isempty(a(l).label)
        total_label_found = total_label_found + 1;
    end
end

if total_label_found == 0
    fprintf('\n> No labels found \n');
else
    fprintf('\n> Labels found in %d atlas(es):\n',total_label_found);
end

for l=1:length(a)
    if ~isempty(a(l).label) 
        fprintf('\n%s',a(l).name)
        for j=1:length(a(l).label)
            fprintf('\n\t%s',a(l).label{j})
        end
        fprintf('\n');
    end
end
        
fprintf('\n');
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
