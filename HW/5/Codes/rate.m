function f = rate(x)
    f = -1 * ( qfunc((x(:,1)-10)/8) + qfunc((x(:,2)-10)/8) + qfunc((x(:,3)-10)/8) + qfunc((x(:,4)-10)/8) );
end
