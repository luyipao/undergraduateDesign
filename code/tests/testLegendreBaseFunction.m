a = -1;
b = 1;

x = [-1 0 1];
disp(['input: ' 'a = ' num2str(a) '; b = ' num2str(b)]);
for i = 0:3
    [PN, DPN] = legendreBaseFunction(i,-1,1);
    disp(['----------' 'n = ' num2str(i) '----------']);
    disp(['x = ' num2str(x)]);
    disp(['PN(x) = ' num2str(PN(x))]);
    disp(['PN(x) = ' num2str(DPN(x))]);
end
