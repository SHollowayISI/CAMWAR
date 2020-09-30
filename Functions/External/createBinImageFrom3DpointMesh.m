function [OutputBinaryImage] = createBinImageFrom3DpointMesh(ShapeCoords, OutputImageSize, TriangularMeshConnectingIndex)
% createBinImageFrom3DpointMesh Computes binary image from 3D mesh-pointset
%   It calculates positions of white pixel within each triangle according
%   to given output image size. Then fills that position with white pixel.
%   After that, it fills hole using inbuilt matlab function 'imfill'.
%   
%   Author: Saurabh J. Shigwan
%   Contact: saurabh.shigwan@cse.iitb.ac.in
%
%   Inputs: 
%   ShapeCoords                   --> (3 x N) matrix, where N is 
%                                     number of points in Mesh.
%
%   TriangularMeshConnectingIndex --> (M x 3) matrix, where M is
%                                     number of triangles in 
%
%   OutputImageSize               --> (1 x 3) row vector, indicating size
%                                     of outputImage
%   
%   Outputs:
%   OutputBinaryImage             --> output binary image of 
%                                     'OutputImageSize' size
%   
%  EXAMPLE is provided in 'demoTest.m' with 'test.mat' file.
%
    ShapeCoords = ShapeCoords([2 1 3],:);
    X0 = ShapeCoords;
    [noOfDims,~] = size(X0);
    t = TriangularMeshConnectingIndex';
    noOfTriangles = size(t,2);
    
    X1 = X0(:,t);
    X1 = reshape(X1,noOfDims,3,noOfTriangles);
    
    size_mat1 = OutputImageSize(1); size_mat2 = OutputImageSize(2); size_mat3 = OutputImageSize(3);
    %size_mat = 1;
    OutputBinaryImage = false(size_mat1, size_mat2, size_mat3);
    pixelWid = 0.3;

    V1 = X1(:,[2:end,1],:) - X1;
    % Select singular triangles 
    singularTriIdx = logical(squeeze(sum(sum(V1.^2,1) == 0,2)));
    
    % Remove Singular triangles
    X1 = X1(:,:,~singularTriIdx);
    V1 = X1(:,[2:end,1],:) - X1;
    
    %Get correct ordering of points
    while true
        t1 = squeeze(-sum(V1(:,1,:).*V1(:,3,:),1) ./ sum( V1(:,1,:).*V1(:,1,:),1));
        t1ValidIdx = (t1 >= 0) & (t1 <= 1);
        
        if sum(~t1ValidIdx) == 0
            break;
        end
        X1(:,:,~t1ValidIdx) = X1(:,[2:end,1],~t1ValidIdx);
        V1 = X1(:,[2:end,1],:) - X1;
    end
    
    C11 = -sum(V1(:,1,:) .* V1(:,1,:),1) ./ sum( V1(:,3,:) .* V1(:,1,:),1);
    C21 = sum( V1(:,1,:) .* V1(:,1,:),1) ./ sum( V1(:,2,:) .* V1(:,1,:),1);
        
    % ========================= First Half Triangle =====================================
    step1 = squeeze(min( pixelWid ./ abs(V1(:,1,:)),[],1));
    maxNoOfsteps1 = max(floor(t1./step1));
    stepTemp1 = t1/maxNoOfsteps1;
    
    t1_list = reshape(bsxfun(@times,1:maxNoOfsteps1,stepTemp1)',1,maxNoOfsteps1,noOfTriangles);
    ptsList = cat(2,X1(:,1,:), X1(:,2,:));
    ptsList = reshape(ptsList,noOfDims,2*noOfTriangles);
    
    t2 = bsxfun(@times,t1_list,C11);
    
    P1 = bsxfun(@plus , X1(:,1,:) , bsxfun( @times,t1_list,V1(:,1,:) ) );
    P2 = bsxfun(@plus , X1(:,1,:) , bsxfun(@times,t2,(-V1(:,3,:) ) ) );
    VP1 = P2-P1;
    stepInner1 = min(min(min((pixelWid)./abs(VP1),[],1),[],2));
    t21 = 0:stepInner1:1;
    t21Length = length(t21);
    VP11 = reshape(VP1,noOfDims,1,maxNoOfsteps1,noOfTriangles);
    P11 = reshape(P1,noOfDims,1,maxNoOfsteps1,noOfTriangles);
    P_list1 = bsxfun(@plus,P11,bsxfun(@times,t21,VP11 ));
    
    ptsList = [ptsList,reshape(P_list1,noOfDims,t21Length*maxNoOfsteps1*noOfTriangles)];
    
    % ========================= Second Half Triangle =====================================
    maxNoOfsteps2 = max(floor((1-t1)./step1));
    stepTemp2 = (1-t1)/maxNoOfsteps2;
    t3_list = reshape(bsxfun(@plus,t1,bsxfun(@times,1:(maxNoOfsteps2-1),stepTemp2))',1,maxNoOfsteps2-1,noOfTriangles);
    
    t4 = bsxfun(@times,-(1-t3_list),C21);
    P3 = bsxfun(@plus , X1(:,1,:) , bsxfun( @times,t3_list,V1(:,1,:) ) );
    P4 = bsxfun(@plus , X1(:,2,:) , bsxfun(@times,t4,V1(:,2,:) ) );
    VP3 = P4-P3;
    stepInner2 = min(min(min((pixelWid)./abs(VP3),[],1),[],2));
    t41 = 0:stepInner2:1;
    t41Length = length(t41);
    VP31 = reshape(VP3,noOfDims,1,maxNoOfsteps2-1,noOfTriangles);
    P31 = reshape(P3,noOfDims,1,maxNoOfsteps2-1,noOfTriangles);
    P_list2 = bsxfun(@plus,P31,bsxfun(@times,t41,VP31 ));
    
    ptsList = [ptsList,reshape(P_list2,noOfDims,t41Length*(maxNoOfsteps2-1)*noOfTriangles)];
    
    ptsList = max( 1, bsxfun(@min,[size_mat1;size_mat2;size_mat3], round(ptsList) ) );
    ptsVec = ptsList(1, :) + (ptsList(2, :) - 1)*size_mat1 + (ptsList(3, :) - 1)*size_mat1*size_mat2;
    
    OutputBinaryImage(ptsVec) = 1; % Fill position with White Pixel
    OutputBinaryImage = imfill(OutputBinaryImage,'holes'); % Fill hole
end