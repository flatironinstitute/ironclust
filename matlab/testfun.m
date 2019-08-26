function a=testfun()

a=test_()

end


function a=test_()
% a=dbstack(1); %, '-completenames');
a=dbstack(1); a=a.name;
end