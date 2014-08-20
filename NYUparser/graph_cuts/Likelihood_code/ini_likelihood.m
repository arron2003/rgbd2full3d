function [ trimap ] = ini_likelihood( user_input,settings )

if(settings.fixed_foreground == 0)
    trimap.B = find(user_input <= 64);
    trimap.U = find(user_input >= 128);
else
    trimap.B = find(user_input <= 64);
    trimap.U = find(user_input == 128);
    trimap.F = find(user_input > 128);
end


end
