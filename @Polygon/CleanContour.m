function contour = CleanContour(PG,cont_original,basepos)
contour = cell(0);
for i =1:numel(cont_original)
    s1 = cont_original{i}(1,1:1:end);
    s2 = cont_original{i}(2,1:1:end);

    cat = 0;
    for j = 1:length(s1)
        [fs1,fs2] = PG.get('2Pos',s1(j),s2(j));
        try
            s2_new = correct_position(PG,s1(j),s2(j),sig);
            if s2_new>PG.S(PG.VL(end))||s2_new<0
                s1_new = correct_position(PG,s2(j),s1(j),sig);
                if s1_new>PG.S(PG.VL(end))||s1_new<0
                    cat = cat+1;
                else
                    s1(j) = s1_new;
                end
            else
                s2(j) = s2_new;
            end
        catch
        end
    end
    theta = zeros(1,size(s1,2));
    d = zeros(2,size(s1,2));
    h = zeros(1,size(s1,2));
    for j = 1:size(s1,2)
        [h(1,j),trans] = PG.get('comLocationDS',basepos,s1(j),s2(j));
        theta(1,j) = trans{1};
        d(:,j) = trans{2};
    end
    if size(s1,2)>3
        contour{numel(contour)+1,1} = [s1;s2;h;theta;d];
    end
end
end