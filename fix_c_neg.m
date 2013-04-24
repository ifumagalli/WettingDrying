function c = fix_c_neg(vertices,elements,c,ar)

    neg_vert = find(c<0);
 
    iv1 = [];
    iv2 = [];
    iv2bis = [];
    iv3 = [];
    for iii = 1:size(elements,2)
        idxs = [find(neg_vert==elements(1,iii)) find(neg_vert==elements(2,iii)) find(neg_vert==elements(3,iii))];
        iv1 = setdiff(negvert(idxs),elements(1:3,iii)); % contains iv1 and iv2 (if iv2 exists)
        if length(iv1)>1
            iv2 = iv1(2);
            iv1 = iv1(1);
            iv3 = negvert(idxs);
        else
            iv3 = negvert(idxs(1));
            iv2bis = negvert(idxs(2));
        end
    
        v1 = vertices(1:2,iv1);
        v2 = vertices(1:2,iv2); % empty if iv2 is empty, i.e. only one vertice is wet
        v2bis = vertices(1:2,iv2bis); % empty if iv2bis is empty, i.e. two vertices are wet
        v3= vertices(1:2,iv3);
        
        xy = [vertices(1,iv1)  vertices(2,iv1)  1;
              vertices(1,iv2)  vertices(2,iv2)  1;
              vertices(1,iv3)  vertices(2,iv3)  1];
        coeff_p = xy\[c(iv1); c(iv2); c(iv3)]; % coeff(1)*x+coeff(2)*y+coeff(3) = c(x,y)
        
        coeff_z1 = [vertices(1,iv1) vertices(2,iv1);
                    vertices(1,iv3) vertices(2,iv3)]   \   [1; 1];
        z1 = [coeff_p(1:2)'; coeff_z1']\[-coeff_p(3);1];

        if isempty(v2)
            coeff_z2 = [vertices(1,iv1) vertices(2,iv1);
                        vertices(1,iv2) vertices(2,iv2)]   \   [1; 1];
            z2 = [coeff_p(1:2)'; coeff_z2']\[-coeff_p(3);1];
            
            aaa = [z1-v1 0; z2-v1 0; 0 0 c(iv1)];
            vol = det(aaa)/6;
        else
            coeff_z2 = [vertices(1,iv2) vertices(2,iv2);
                    vertices(1,iv3) vertices(2,iv3)]   \   [1; 1];
            z2 = [coeff_p(1:2)'; coeff_z2']\[-coeff_p(3);1];
        
            aaa = [z1-v1 0; z2-v1 0; 0 0 c(iv1)];
            bbb = [z1-v2 0; v1-v2 0; 0 0 c(iv2)];
            w1 = [0 0 c(iv2)] - v1;
            a1 = norm(z1-v1,2); a2 = norm(z1-v2); a3 = norm(v1-v2);
            p = (a1+a2+a3)/2;
            area = sqrt(p*(p-a1)*(p-a2)*(p-a3));
            altezza = area/norm(v1-v2);
            proiez = sqrt(norm(v2-z1)^2 - altezza^2);
            w2 = [z1-proiez*(v1-v2)/norm(v1-v2) 0];
            ccc = [0 0 c(iv1); w1; w2];
            vol = (det(aaa)+det(bbb)+det(ccc))/6;
        end
        
        
    end
    
    % rescaling heights
    c(iv3) = 0;
    if isempty(v2)
        c(iv2bis) = 0;
        area = pdetrg(vertices(1:2,iii),elements(1:3,iii));
        c(iv1) = 6*vol/area;
    else
        c(iv1) = vol/((norm(v2-v1))^2);
        c(iv2) = c(iv1);
    end
end