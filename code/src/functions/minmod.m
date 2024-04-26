
function y = minmod(v,M)
	index = abs(v(:,1)) <= M;
	y(index) = v(index,1);
	y(~index) = all(sign(v(~index,:)) == sign(v(~index,1)),2) .* sign(v(~index,1)) .* min(abs(v(~index)),[],2);
end