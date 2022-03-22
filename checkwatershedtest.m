function checkWater = checkwatershedtest(t1,t2,tw)
checkWater = false;
% this function checks if there is a high watershed curve between A and B
% this is only relevant after translation, and before checking rotation
% s = A(1);
% contact_pos = PG.get('1pos',s);
% contact_to_com = PG.COM-contact_pos;
% thetaWater = atan2(contact_to_com(1),contact_to_com(2));

Delta1 = wrapToPi(t2-t1);
Delta2 = wrapToPi(tw-t1);
if Delta1*Delta2>0 && abs(Delta1)>abs(Delta2)
    checkWater = true;
end
end