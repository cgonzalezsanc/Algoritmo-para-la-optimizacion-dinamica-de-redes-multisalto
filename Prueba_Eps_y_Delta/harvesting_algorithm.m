function [e_op, e_out, e_in] = harvesting_algorithm(p, h, c_p, B_cap, B_state, pi)
% Algorithm that calculates the power taken out from the battery, bought
% from the operator and harvested. In this case there is no cost associated
% to introduce energy in or take it out from the battery.
%
% p = Power required
% h = Power harvested
% c_p = Power cost
% B_cap = Battery cap
% B_state = Battery state
% pi = lagrangian multiplier

% Case where getting the power from the operator is the best choice
if c_p < pi
    e_in = min(h, B_cap-B_state); %check if we exceed the maximum capacity of the battery. If so, we use the energy overflowed
    e_out = 0;
    e_op = max(0,p-(h-e_in));
    return
end

% We get the power from the battery.
% Case where we can provide enough energy from harvesting
if p <= h
    e_in = min(h-p, B_cap-B_state);
    e_out = 0;
    e_op = 0;
    return
end

% Case when we also need to take energy out from the battery
e_in = 0;
e_out = min(p-h, B_state); % check if we have enough energy
e_op = p-h-e_out; % we still may need to buy energy from the operator

end

