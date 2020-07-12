%2 Dimensional Histogram
%Counts number of points in the bins defined by vYEdge, vXEdge.
%size(vX) == size(vY) == [n,1]
%size(mHist) == [length(vYEdge) -1, length(vXEdge) -1]
%
%EXAMPLE
%   mYX = rand(100,2);
%   vXEdge = linspace(0,1,10);
%   vYEdge = linspace(0,1,20);
%   mHist2d = hist2d(mYX,vYEdge,vXEdge);
%
%   nXBins = length(vXEdge);
%   nYBins = length(vYEdge);
%   vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins));
%   vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins));
%   pcolor(vXLabel, vYLabel,mHist2d); colorbar
function [mHist, Index] = hist2d (Vel, Rho, vYEdge, vXEdge)
%Rho=reshape(RhoMdlI,1,length(YVel)*length(ZVel))';
%Vel=reshape(VelMdl,1,length(YVel)*length(ZVel))';

%nRho=length(Rho)-sum(isnan(Rho));

%Define Rho (X) and Vel (Y) bins
%vXEdge = linspace(0.5,3.5,16);  
%vYEdge = linspace(4.5,8.5,21);  


nRow = length (vYEdge)-1;
nCol = length (vXEdge)-1;

%Initialize matrix to store occurance
mHist = zeros(nRow,nCol);
RowIndex = zeros(1,length(Rho));
%Index = zeros(1,2);

for iRow = 1:nRow
    Index{iRow} = zeros(1,2);
    rRowLB = vYEdge(iRow);
    rRowUB = vYEdge(iRow+1);
 
    %Rho values correlated to vel values within the first row of bins
    RhoFound = Rho((Vel > rRowLB) & (Vel <= rRowUB)); % look for Vel values between lower and upper bounds
    RowIndex(1,:) = ((Vel > rRowLB) & (Vel <= rRowUB));
    RhoIndex = find(RowIndex~=0)'; % find which rho's that these velocities correspond with      
    if (~isempty(RhoFound))       
        %Counts the Rho values found in column bins of the first row
        [vFound Bin] = histc (RhoFound, vXEdge);       %vFound tells number of rhos within each vel bin 
        
        %Store the occurrance values in mHist matrix
        mHist(iRow, :)= vFound(1:length(vFound)-1)'; % store number of occurrences in matrix
        %Store the Bin index and RhoIndex in Index
        temp=[Bin+(iRow-1)*nCol RhoIndex];
        Index{iRow}=cat(1,Index{iRow},temp);
        clear temp;
        %for j=1:length(RhoIndex)
        %    Index(min(find(Index==0))+j-1,1)=Bin(j);
        %    Index(min(find(Index==0))+j-1,2)=RhoIndex(j);
        %end
    end
    Index{iRow}(1,:)=[];
end
%Index(1,:)=[];

%Plot Histogram
% nXBins = length(vXEdge);  
% nYBins = length(vYEdge);  
% vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins));  
% vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins));  
% pcolor(vXLabel, vYLabel,mHist); colorbar
% colormap(flipud(bone));
