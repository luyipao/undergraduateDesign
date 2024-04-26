% testDividedDiff
x = [0 1 2 3];
y = [6 -3 -6 9];
x3 =  x1(1:2);
y3 = y1(1:2);
table = dividedDiff(x,y);
table3 = dividedDiff(x3,y3);
table2 = dividedDiff([2 3],[-6 9],table3);
e = table - table2;