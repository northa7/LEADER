% Phase correction, if the phase angle changes too much
% (limits: crit_change1 if the phase angle < crit_low,
% otherwise crit_change2)
crit_change1 = 1;
crit_change2 = 2;
crit_low = 8;
for i=1:Nappar
    % The start and end indices of the data
    ind = sum(pointsperapp(1:i-1))+1;
    inde = ind+pointsperapp(i)-1;
    temp = ang_back(ind:inde);
    
    % Do we need exponential or linear correction
    expcorr = ( min(temp) < crit_low ) && ...
        ( max(temp) - min(temp) > crit_change1 );
    lincorr = ( min(temp) >= crit_low ) && ...
        ( max(temp) - min(temp) > crit_change2 );
    
    if ( expcorr )
        % Phase correction with an exponential fit
        AA = [temp.^0, temp];
        yy = log(L_back(ind:inde));
        xx = AA\yy;
        % Debug: make sure we have a decreasing exponential curve
        if ( xx(2) > 0 )
            XX = [1; 0];
        else
            XX = [exp(xx(1)); xx(2)];
        end
        L_back(ind:inde) = L_back(ind:inde) ./ ...
            (XX(1)*exp(XX(2)*temp));
        clear XX
            
    elseif ( lincorr )
        % Phase correction with a linear curve
        AA = [temp, temp.^0];
        yy = L_back(ind:inde);
        xx = AA\yy;
        % Debug: make sure we have a decreasing linear curve
        if ( xx(1) > 0 )
            XX = [0; 1];
        else
            XX = xx;
        end
        L_back(ind:inde) = L_back(ind:inde) ./ ...
            (XX(1)*temp+XX(2));
        clear XX
    end
end

% If we didn't get enough points, combine the apparitions
if ( Nappar > 1 && any(pointsperapp < wanted) )
    % Set the scaling according to the 1st apparition
    avg_scale = mean(L_back(1:pointsperapp(1)));
    % Scale the other apparitions to this
    for i=2:Nappar
        % The start and end indices of the data
        ind = sum(pointsperapp(1:i-1))+1;
        inde = ind+pointsperapp(i)-1;
        avg_temp = mean(L_back(ind:inde));
        scale_factor = avg_scale/avg_temp;
        % Re-scale the apparition
        L_back(ind:inde) = scale_factor * L_back(ind:inde);
    end
    % After combining, we have only 1 effective apparition
    Nappar_eff = 1;
end