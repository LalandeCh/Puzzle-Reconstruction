function Main

%% Algorithm of puzzle reconstruction
%% Optimization Module

%% Lalande Chatain Benjamin
%% Katherine Sheran
%% Carmen Moreno Genis



%% Initialization & Download of the Data of each Fragments

clear all, close all, clc

Puzzle1 = load('HK-233_a_09-Sep-2015_ver1_XYZ+_newZ_nb_clean.mat');

Puzzle2 = load('HK-233_b_09-Sep-2015_ver1_XYZ+_newZ_nb_clean.mat');

Puzzle3 = load('HK-233_c_09-Sep-2015_ver1_XYZ+_newZ_nb_clean.mat');

%% Transformation of the Data to 3D Point Clouds

ptCloud1 = Create_PtCloud (Puzzle1.Xclean(:), Puzzle1.Yclean(:), Puzzle1.Zclean(:));
ptCloud2 = Create_PtCloud (Puzzle2.Xclean(:), Puzzle2.Yclean(:), Puzzle2.Zclean(:));
ptCloud3 = Create_PtCloud (Puzzle3.Xclean(:), Puzzle3.Yclean(:), Puzzle3.Zclean(:));

%% Method of Facets Segmentation Based on Matlab Functions

%   Sub_PtClouds = Facet_Segmentation_PCFitPlane (ptCloud3, 2);
%   ShowSubPtCloud (Sub_PtClouds)

%% Proposed Method of Facets Segmentation

% Sub_PtClouds3 = Facet_Segmentation_Proposed (ptCloud3, 1, 10);

%% Save Facets Segmented of each Fragments
%% Do Not Uncomment

% save('Sub_PtClouds1.mat', 'Sub_PtClouds1') 
% save('Sub_PtClouds2.mat', 'Sub_PtClouds2') 
% save('Sub_PtClouds3.mat', 'Sub_PtClouds3')

%% Load and display each Facets of a Fragment and the Reconstruct Fragment from the Facets 

load('Sub_PtClouds1')
% ShowSubPtCloud (Sub_PtClouds1)

load('Sub_PtClouds2')
% ShowSubPtCloud (Sub_PtClouds2)

load('Sub_PtClouds3')
% ShowSubPtCloud (Sub_PtClouds3)

%% Find Matches between the Facets of 2 Fragments and perfom a merging

[Ind1, Ind2, tform_Matrix] = FindMatchFacets(Sub_PtClouds2, Sub_PtClouds3);

MergingPtClouds (ptCloud2, ptCloud3, Ind1, Ind2, tform_Matrix);


end

function ShowSubPtCloud (Sub_PtClouds)

% This function shows the different facets of a fragment separately, then
% the reconstructed fragment using a color for each facet in order to see
% the result of the segmentation

DataBase_Color = [0 0 0; 255 0 0; 255 255 0; 153 255 51; 0 102 0; 0 255 255; 0 0 255; 153 0 153; 255 0 127];
figure

for i = 1 : size(Sub_PtClouds, 2)
    
    [X, Y] = size (Sub_PtClouds{i}.Location);
    
    CodeColor = uint8(zeros(X, Y));
    CodeColor(:,1)=DataBase_Color(i, 1);
    CodeColor(:,2)=DataBase_Color(i, 2);
    CodeColor(:,3)=DataBase_Color(i, 3);
    
    Sub_PtClouds{i} = pointCloud(Sub_PtClouds{i}.Location, 'Color', CodeColor);
    
    subplot(1, size(Sub_PtClouds, 2), i)
    
    pcshow(Sub_PtClouds{i})
    
    Sub_PtClouds{1} = pcmerge(Sub_PtClouds{1}, Sub_PtClouds{i}, 1);
    
end

figure
pcshow(Sub_PtClouds{1})

end

function ptCloud = Create_PtCloud ( XClean, YClean, ZClean )

Convex_Alpha_Shape = alphaShape( XClean, YClean, ZClean, 5);

ptCloud = pointCloud( Convex_Alpha_Shape.Points );

end

function Sub_PtClouds = Facet_Segmentation_Proposed(ptCloud, Threshold, PercentageFacet)

% This function perform the segmentation of a fragment into several facets.

Sub_PtClouds = [];
inlierIndices= [];
outlierIndices= [];

% RemainPtCloud is a copy of the input Point Could on which we will perform
% the method

RemainPtCloud = ptCloud;
Ind = 1;

disp('Initialization of the facet segmentation')

% Definition of the rate of lost of information to 10 %

while (RemainPtCloud.Count > ptCloud.Count / 10)
    
    disp('New search of facet')
    IndWhile = 0;
    
    inlierIndices= [];
    outlierIndices= [];
    
    % Definition of the minimum numbers of points of a facet comparing to
    % the complete original Point Cloud
    
    while ((size(inlierIndices, 2) < ptCloud.Count * (PercentageFacet / 100)))
        
        inlierIndices= [];
        outlierIndices= [];
        
        % Select the index of 2 random points from the point cloud
        
        IndP(1) = round(rand() * RemainPtCloud.Count);
        IndP(2) = round(rand() * RemainPtCloud.Count);
        
        % On the plane xy, find the coefficients of the equation of the line
        % corresponding to these points
        
        CoeffL = polyfit(RemainPtCloud.Location(IndP, 1), RemainPtCloud.Location(IndP, 2), 1);
        
        % Select 2 points from the line equation
        
        PointLine1 = [-10, CoeffL(1) * -10 + CoeffL(2)];
        PointLine2 = [10, CoeffL(1) * 10 + CoeffL(2)];
        
        A = PointLine1 - PointLine2;
        
        disp('Compute Distance Points to line')
        
        % Compute for each point of the PointCloud the Orthogonal Distance
        % to the line in order to determine if this line can characterise
        % a facet
        
        for i = 1 : RemainPtCloud.Count
            
            B = [RemainPtCloud.Location(i, 1), RemainPtCloud.Location(i, 2)] - PointLine2;
            OrthoD = norm(det([A; B])) / norm(A);
            
            % If the orthogonal distance is inferior to the treshold, we
            % consider the point as inlier
            
            if (OrthoD < Threshold)
                inlierIndices = [inlierIndices i];
            else
                outlierIndices = [outlierIndices i];
            end
            
        end
       
        IndWhile = IndWhile +1;
        
        % If the algorthm is computing a new line more than 50 times
        % without finfing a results, the algrithm stop
        
        if (IndWhile > 50)
            break
        end
        
    end
    
    if (IndWhile > 50)
            break
    end
        
    disp('Finding Facet')
    
    % Select all the inliers to create a Point Cloud representing a facet
    % The outliers form a Point Cloud in which we will attempt to extract 
    % an another facet
    
    Sub_PtClouds{Ind} = select(RemainPtCloud, inlierIndices);
    PT = select(RemainPtCloud, outlierIndices);
    RemainPtCloud = PT;
    
    disp('Count RemainPtCloud')
    display(RemainPtCloud.Count)

    Ind = Ind +1;
    
end

end

function Sub_PtClouds = Facet_Segmentation_PCFitPlane (ptCloud, MaxDistance)

% This function is based on an example presented on the website MathWorks
% https://fr.mathworks.com/help/vision/ref/pcfitplane.html

Pt_Count = ptCloud.Count;
RemainPtCloud_Count = ptCloud.Count;
RemainPtCloud = ptCloud;
Ind = 1;

Sub_PtClouds={};

while (RemainPtCloud_Count > (Pt_Count / 20))
    
    % pcfitplane define a plane on the image, then determine if the points 
    % are inlier or outlier
    
    [model , inlierIndices,outlierIndices] = pcfitplane(RemainPtCloud, 3);
    
    % Select all the inliers to create a Point Cloud representing a facet
    % The outliers form a Point Cloud in which we will attempt to extract 
    % an another facet
    
    Sub_PtClouds{Ind} = select(RemainPtCloud, inlierIndices);
    RemainPtCloud = select(RemainPtCloud, outlierIndices);
    
    RemainPtCloud_Count = RemainPtCloud.Count;
        
    Ind = Ind + 1
    
end

end

function  [Ind_Sub_PtClouds2,Ind_Sub_PtClouds1, tform_Matrix] = FindMatchFacets(Sub_PtClouds1, Sub_PtClouds2)

for i = 1 : length(Sub_PtClouds1)    
  for j = 1 :  length(Sub_PtClouds2)
    
      % Considering 2 fragments, compare their facets between each others
      % using the function pcregrigid. The output is the root mean square
      % error and the matrix transformation between the 2 facets
      
          [tform, ~, rmse] = pcregrigid(Sub_PtClouds1{i}, Sub_PtClouds2{j});
          
          tform_Matrix{i, j} = tform;
          rmse_Matrix(i, j) = rmse;
        
  end
end

% Extrat the minimum value of the rmse corresponding to best matches
% between facets

[Ind_Sub_PtClouds1, Ind_Sub_PtClouds2] = find(rmse_Matrix == min(min(rmse_Matrix)));

end

function MergingPtClouds (ptCloud1, ptCloud2, Ind1, Ind2, tform_Matrix)

% Apply the matrix transformation corresponding to the different Point
% Clouds, then merge the point clouds.

ptCloud1_Transform = pctransform(ptCloud1, tform_Matrix{Ind1, Ind2});

figure
pcshow(pcmerge(ptCloud1_Transform, ptCloud2, 1));

end



