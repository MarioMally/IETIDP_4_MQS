function [nrbarrCond2D, nrbarrAir2D, nrbarrCond3D, nrbarrAir3D] = planar_coil_mario(thicknessCond, thicknessAir, midOffset ,cornerRad, numQuarterTurn, airBox)

    % Rotation matrix
    R = @(alpha) [cosd(alpha),-sind(alpha);...
                  sind(alpha),cosd(alpha)];
    
    %% Parameter setup
    innerRad = cornerRad;
    outerRad = thicknessCond + cornerRad;
    innerVec = [innerRad; 0];
    outerVec = [outerRad; 0];
    countAir = 0;
    countCond = 0;
    alpha = 0;
    
    %% Coil interior start
    p1(:,1) = [thicknessCond/2;0]; % Start of coil side 1
    p1(:,end+1) = p1(:,end) + [thicknessCond; 0];
    p1(:,end+1) = p1(:,end) + [midOffset - thicknessCond - innerRad - thicknessCond/2; 0];
    
    p2(:,1) = p1(:,1) + [0; thicknessCond];
    p2(:,end+1) = p2(:,end) + [thicknessCond; 0];
    p2(:,end+1) = p2(:,end) + [midOffset - thicknessCond - innerRad - thicknessCond/2; 0];
    
    center = p1(:,end) + R(alpha-90)*innerVec;
    p1(:,end+1) = center + R(alpha+45)*innerVec;
    p2(:,end+1) = center + R(alpha+45)*outerVec;
    
    nrbarrCond2D(1) = nrb4surf(p1(:,1), p1(:,2), p2(:,1), p2(:,2));
    nrbarrCond2D(2) = nrb4surf(p1(:,2), p1(:,3), p2(:,2), p2(:,3));
    
    innerArc = nrbcirc(innerRad, center, deg2rad(alpha+45), deg2rad(alpha+90));
    outerArc = nrbcirc(outerRad, center, deg2rad(alpha+45), deg2rad(alpha+90));
    nrbarrCond2D(end+1) = nrbreverse(nrbruled(innerArc,outerArc),1);
    
    %% Quarter turns of coil
    for i=1:numQuarterTurn
        
        p1(:,end+1) = center + R(alpha)*innerVec;
        p2(:,end+1) = center + R(alpha)*outerVec;
    
        innerArc = nrbcirc(innerRad, center, deg2rad(alpha+0), deg2rad(alpha+45));
        outerArc = nrbcirc(outerRad, center, deg2rad(alpha+0), deg2rad(alpha+45));
        nrbarrCond2D(end+1) = nrbreverse(nrbruled(innerArc,outerArc),1);
        
        lenVec1 = [midOffset + countAir*thicknessAir + countCond*thicknessCond; 0];
        lenVec2 = [midOffset + countAir*thicknessAir + (countCond+1)*thicknessCond; 0];
        p1(:,end+1) = p1(:,end) + R(alpha-90)*lenVec1 - 2*R(alpha-90)*innerVec;
        p2(:,end+1) = p2(:,end) + R(alpha-90)*lenVec2 - R(alpha-90)*innerVec - R(alpha-90)*outerVec;
    
        nrbarrCond2D(end+1) = nrb4surf(p1(:,end-1), p1(:,end), p2(:,end-1), p2(:,end));
    
        alpha = alpha - 90;
        center = p1(:,end) + R(alpha-90)*innerVec;
        p1(:,end+1) = center + R(alpha+45)*innerVec;
        p2(:,end+1) = center + R(alpha+45)*outerVec;
    
        innerArc = nrbcirc(innerRad, center, deg2rad(alpha+45), deg2rad(alpha+90));
        outerArc = nrbcirc(outerRad, center, deg2rad(alpha+45), deg2rad(alpha+90));
        nrbarrCond2D(end+1) = nrbreverse(nrbruled(innerArc,outerArc),1);
        
        if rem(i,2)==1
            countAir = countAir + 1;
        else
            countCond = countCond + 1;
        end
    
    end
    
    %% Endpart coil
    innerArc = nrbcirc(innerRad, center, deg2rad(alpha+0), deg2rad(alpha+45));
    outerArc = nrbcirc(outerRad, center, deg2rad(alpha+0), deg2rad(alpha+45));
    nrbarrCond2D(end+1) = nrbreverse(nrbruled(innerArc,outerArc),1);
    
    p1(:,end+1) = center + R(alpha)*innerVec;
    p2(:,end+1) = center + R(alpha)*outerVec;
    lenVec1 = [midOffset + countAir*thicknessAir + countCond*thicknessCond - 6*thicknessCond; 0];
    lenVec2 = [midOffset + countAir*thicknessAir + (countCond+1)*thicknessCond - 6*thicknessCond; 0];
    p1(:,end+1) = p1(:,end) + R(alpha-90)*[2;0] - 2*R(alpha-90)*innerVec;
    p2(:,end+1) = p2(:,end) + R(alpha-90)*[2;0] - R(alpha-90)*innerVec - R(alpha-90)*outerVec;
    
    nrbarrCond2D(end+1) = nrb4surf(p1(:,end-1), p1(:,end), p2(:,end-1), p2(:,end));
    
    lenVec = [thicknessCond; 0];
    p1(:,end+1) = p1(:,end) + R(alpha-90)*lenVec;
    p2(:,end+1) = p2(:,end) + R(alpha-90)*lenVec;
    nrbarrCond2D(end+1) = nrb4surf(p1(:,end-1), p1(:,end), p2(:,end-1), p2(:,end));
    
    % Shift midpoint of coil
    p1 = p1 + repmat([-midOffset/2;midOffset/2],1,size(p1,2));
    p2 = p2 + repmat([-midOffset/2;midOffset/2],1,size(p2,2));
    for i=1:numel(nrbarrCond2D)
        nrbarrCond2D(i) = nrbtform(nrbarrCond2D(i), vectrans([-midOffset/2,midOffset/2,0]));
    end
    
    %% Air at interior corner
    nrbarrAir2D(1) = nrbruled(nrbextract(nrbarrCond2D(12),3), nrbextract(nrbarrCond2D(1),1));
    nrbarrAir2D(2) = nrbruled(nrbextract(nrbarrCond2D(13),3), nrbextract(nrbarrCond2D(1),4));
    
    %% Air inside of coil
    p13 = [-midOffset; -midOffset]/3 - [thicknessCond; 0];
    p14 = [midOffset; -midOffset]/3;
    p15 = [-midOffset; midOffset]/3 - [thicknessCond; 0];
    p16 = [midOffset; midOffset]/3;
    nrbarrAir2D(3) = nrb4surf(p13,p14,p15,p16);
    
    nrbarrAir2D(4) = nrbruled(nrbextract(nrbarrAir2D(3), 4), nrbextract(nrbarrCond2D(2), 3));
    nrbarrAir2D(5) = nrbruled(nrbreverse(nrbextract(nrbarrAir2D(3), 2), 1), nrbextract(nrbarrCond2D(5), 3));
    nrbarrAir2D(6) = nrbruled(nrbreverse(nrbextract(nrbarrAir2D(3), 3), 1), nrbextract(nrbarrCond2D(8), 3));
    nrbarrAir2D(7) = nrbruled(nrbextract(nrbarrAir2D(3), 1), nrbextract(nrbarrCond2D(11), 3));
    
    nrbarrAir2D(8) = nrbcoons(nrbextract(nrbarrAir2D(5), 1),...
                            nrbextract(nrbarrCond2D(3), 3),...
                            nrbextract(nrbarrAir2D(4), 2),...
                            nrbreverse(nrbextract(nrbarrCond2D(4), 3), 1));
    nrbarrAir2D(9) = nrbcoons(nrbextract(nrbarrAir2D(5), 2),...
                            nrbreverse(nrbextract(nrbarrCond2D(7), 3), 1),...
                            nrbextract(nrbarrAir2D(6), 1),...
                            nrbextract(nrbarrCond2D(6), 3));
    nrbarrAir2D(10) = nrbcoons(nrbextract(nrbarrAir2D(6), 2),...
                            nrbreverse(nrbextract(nrbarrCond2D(10), 3), 1),...
                            nrbextract(nrbarrAir2D(7), 1),...
                            nrbextract(nrbarrCond2D(9), 3));
    nrbarrAir2D(11) = nrbcoons(nrbextract(nrbarrAir2D(7), 2),...
                            nrbreverse(nrbextract(nrbarrCond2D(1), 3), 1),...
                            nrbextract(nrbarrAir2D(4), 1),...
                            nrbextract(nrbarrAir2D(1), 1));
    
    %% Air in between conductor strands
    for i=1:numel(nrbarrCond2D)-13
        nrbarrAir2D(end+1) = nrbruled(nrbextract(nrbarrCond2D(i+13),3), nrbextract(nrbarrCond2D(i+1),4));
    end
    
    %% Air outside of coil
    box(:,1) = [-airBox/2; -airBox/2];
    box(:,2) = [airBox/2; -airBox/2];
    box(:,3) = [airBox/2; airBox/2];
    box(:,4) = [-airBox/2; airBox/2];
    
    dist2boxCorners = vecnorm(box - repmat(p1(:,end), 1, size(box, 2)));
    [~,shift] = min(dist2boxCorners);
    box = circshift(box, size(box, 2) - shift + 1, 2);
    
    lenVecBox = [airBox; 0];
    
    box(:,5) = box(:,1) + (1/6)*R(alpha+90)*lenVecBox;
    box(:,6) = box(:,1) + (1/3)*R(alpha+90)*lenVecBox;
    box(:,7) = box(:,1) + (2/3)*R(alpha+90)*lenVecBox;
    
    box(:,8) = box(:,2) + (1/3)*R(alpha+180)*lenVecBox;
    box(:,9) = box(:,2) + (2/3)*R(alpha+180)*lenVecBox;
    
    box(:,10) = box(:,3) + (1/3)*R(alpha+270)*lenVecBox;
    box(:,11) = box(:,3) + (2/3)*R(alpha+270)*lenVecBox;
    
    box(:,12) = box(:,4) + (1/3)*R(alpha)*lenVecBox;
    box(:,13) = box(:,4) + (2/3)*R(alpha)*lenVecBox;
    box(:,14) = box(:,4) + (5/6)*R(alpha)*lenVecBox;
    
    box(:,15) = (box(:,14) + box(:,13) + p1(:,end) + p2(:,end-12) + p2(:,end-11) + box(:,1))/6;
    
    nrbarrAir2D(end+1) = nrbcoons(nrbextract(nrbarrCond2D(end-11), 4),...
                                nrbline(p1(:,end), box(:,15)),...
                                nrbreverse(nrbextract(nrbarrAir2D(end), 2), 1),...
                                nrbline(p2(:,end-11), box(:,15)));
    nrbarrAir2D(end+1) = nrbcoons(nrbline(p1(:,end), box(:,1)),...
                                nrbline(box(:,15), box(:,14)),...
                                nrbextract(nrbarrAir2D(end), 4),...
                                nrbline(box(:,1), box(:,14)));
    nrbarrAir2D(end+1) = nrbcoons(nrbextract(nrbarrAir2D(end), 4),...
                                nrbline(p2(:,end-11), box(:,13)),...
                                nrbreverse(nrbextract(nrbarrAir2D(end-1), 2), 1),...
                                nrbline(box(:,14), box(:,13)));
    nrbarrAir2D(end+1) = nrbcoons(nrbextract(nrbarrAir2D(end-1), 3),...
                                nrbline(p2(:,end), box(:,5)),...
                                nrbextract(nrbarrCond2D(end), 2),...
                                nrbline(box(:,1), box(:,5)));
    nrbarrAir2D(end+1) = nrbcoons(nrbextract(nrbarrCond2D(end), 4),...
                                nrbline(box(:,6), box(:,5)),...
                                nrbline(p2(:,end-1), box(:,6)),...
                                nrbextract(nrbarrAir2D(end), 4));
    nrbarrAir2D(end+1) = nrbcoons(nrbextract(nrbarrCond2D(end-1), 4),...
                                nrbline(box(:,7), box(:,6)),...
                                nrbline(p2(:,end-2), box(:,7)),...
                                nrbextract(nrbarrAir2D(end), 1));
    nrbarrAir2D(end+1) = nrbcoons(nrbextract(nrbarrCond2D(end-2), 4),...
                                nrbline(box(:,2), box(:,7)),...
                                nrbline(p2(:,end-3), box(:,2)),...
                                nrbextract(nrbarrAir2D(end), 1));
    nrbarrAir2D(end+1) = nrbcoons(nrbextract(nrbarrCond2D(end-3), 4),...
                                nrbline(box(:,8), box(:,2)),...
                                nrbline(p2(:,end-4), box(:,8)),...
                                nrbextract(nrbarrAir2D(end), 1));
    nrbarrAir2D(end+1) = nrbcoons(nrbextract(nrbarrCond2D(end-4), 4),...
                                nrbline(box(:,9), box(:,8)),...
                                nrbline(p2(:,end-5), box(:,9)),...
                                nrbextract(nrbarrAir2D(end), 1));
    nrbarrAir2D(end+1) = nrbcoons(nrbextract(nrbarrCond2D(end-5), 4),...
                                nrbline(box(:,3), box(:,9)),...
                                nrbline(p2(:,end-6), box(:,3)),...
                                nrbextract(nrbarrAir2D(end), 1));
    nrbarrAir2D(end+1) = nrbcoons(nrbextract(nrbarrCond2D(end-6), 4),...
                                nrbline(box(:,10), box(:,3)),...
                                nrbline(p2(:,end-7), box(:,10)),...
                                nrbextract(nrbarrAir2D(end), 1));
    nrbarrAir2D(end+1) = nrbcoons(nrbextract(nrbarrCond2D(end-7), 4),...
                                nrbline(box(:,11), box(:,10)),...
                                nrbline(p2(:,end-8), box(:,11)),...
                                nrbextract(nrbarrAir2D(end), 1));
    nrbarrAir2D(end+1) = nrbcoons(nrbextract(nrbarrCond2D(end-8), 4),...
                                nrbline(box(:,4), box(:,11)),...
                                nrbline(p2(:,end-9), box(:,4)),...
                                nrbextract(nrbarrAir2D(end), 1));
    nrbarrAir2D(end+1) = nrbcoons(nrbextract(nrbarrCond2D(end-9), 4),...
                                nrbline(box(:,12), box(:,4)),...
                                nrbline(p2(:,end-10), box(:,12)),...
                                nrbextract(nrbarrAir2D(end), 1));
    nrbarrAir2D(end+1) = nrbcoons(nrbextract(nrbarrCond2D(end-10), 4),...
                                nrbline(box(:,13), box(:,12)),...
                                nrbextract(nrbarrAir2D(end-11), 4),...
                                nrbextract(nrbarrAir2D(end), 1));
    %% Equalize Degree of nurbs
    nrbarr2D = [nrbarrCond2D, nrbarrAir2D];
    degarr1 = arrayfun(@(nrb) nrb.order(1), nrbarr2D);
    degarr2 = arrayfun(@(nrb) nrb.order(2), nrbarr2D);
    
    degmax1 = max(degarr1);
    degmax2 = max(degarr2);
    degmax = max([degmax1,degmax2]);
    
    elevarr1 = degmax - degarr1;
    elevarr2 = degmax - degarr2;
    
    for i = 1:numel(nrbarr2D)
        nrbarr2D(i) = nrbdegelev(nrbarr2D(i), [elevarr1(i), elevarr2(i)]);
    end
    
    %% Compute interfaces and boundaries (for checking correctness)
    [interfaces, boundaries] = nrbmultipatch(nrbarr2D);
    
    %% Extend to 3D
    % First extrude coil and air by thickness
    for i=1:numel(nrbarrCond2D)
        nrbarrCond3D(i) = nrbextrude(nrbarrCond2D(i), [0,0,thicknessCond]);
    end
    for i=1:numel(nrbarrAir2D)
        nrbarrAir3D(i) = nrbextrude(nrbarrAir2D(i), [0,0,thicknessCond]);
    end
    % Add layer below coil (with endings)
    for i=1:numel(nrbarrCond2D)
        if i==1
            nrbarrCond3D = [nrbextrude(nrbarrCond2D(i), [0,0,-airBox/2]), nrbarrCond3D];
        elseif i==numel(nrbarrCond2D)
            nrbarrCond3D = [nrbextrude(nrbarrCond2D(i), [0,0,-airBox/2]), nrbarrCond3D];
        else
            nrbarrAir3D(end+1) = nrbextrude(nrbarrCond2D(i), [0,0,-airBox/2]);
        end
    end
    for i=1:numel(nrbarrAir2D)
        nrbarrAir3D(i) = nrbextrude(nrbarrAir2D(i), [0,0,-airBox/2]);
    end
    % Add layer above coil (without endings)
    for i=1:numel(nrbarrCond2D)
        nrbarrAir3D(end+1) = nrbextrude(nrbarrCond2D(i), [0,0,airBox/2 - thicknessCond]);
        nrbarrAir3D(end) = nrbtform(nrbarrAir3D(end), vectrans([0,0,thicknessCond]));
    end
    for i=1:numel(nrbarrAir2D)
        nrbarrAir3D(end+1) = nrbextrude(nrbarrAir2D(i), [0,0,airBox/2 - thicknessCond]);
        nrbarrAir3D(end) = nrbtform(nrbarrAir3D(end), vectrans([0,0,thicknessCond]));
    end
    
    %% Plotting
    % figure(1)
    % clf()
    % plot(p1(1,:),p1(2,:),'Color','b','Marker','x');
    % hold on;
    % plot(p2(1,:),p2(2,:),'Color','r','Marker','x');
    % % plot([-airBox(1)/2,airBox(1)/2,airBox(1)/2,-airBox(1)/2,-airBox(1)/2],...
    % %      [-airBox(2)/2,-airBox(2)/2,airBox(2)/2,airBox(2)/2,-airBox(2)/2],'Color','k')
    % axis equal;
    % grid on;
    
    figure(1)
    clf()
    
    subplot(2,2,1);
    for i=1:numel(nrbarrCond2D)
        mynrbplot(nrbarrCond2D(i),[5,5],'#CC4C03');
        hold on;
    end
    for i=1:numel(nrbarrAir2D)
        mynrbplot(nrbarrAir2D(i),[5,5],'#0080FE');
    end
    axis equal;
    grid on;
    view(0,90);
    xlim([-airBox/2, airBox/2]);
    ylim([-airBox/2, airBox/2]);
    
    subplot(2,2,2);
    for i=1:numel(interfaces)
        patch1 = interfaces(i).patch1;
        side1 = interfaces(i).side1;
        mynrbplot(nrbextract(nrbarr2D(patch1), side1), 3, 'k');
        hold on;
    end
    for i=1:numel(boundaries)
        patch1 = boundaries(i).patches;
        side1 = boundaries(i).faces;
        mynrbplot(nrbextract(nrbarr2D(patch1), side1), 3, 'k');
    end
    axis equal;
    grid on;
    view(0,90);
    xlim([-airBox/2, airBox/2]);
    ylim([-airBox/2, airBox/2]);
    
    subplot(2,2,3)
    for i=1:numel(nrbarrCond3D)
        mynrbplot(nrbarrCond3D(i),[2,0,0],'#CC4C03');
        hold on;
    end
    axis equal;
    xlim([-airBox/2, airBox/2]);
    ylim([-airBox/2, airBox/2]);
    zlim([-airBox/2, airBox/2]);
    view(-207.2911, 44.0191)

end
