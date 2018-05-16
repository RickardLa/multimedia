function output = izigZagTest(input)

% Scans upper  left of matrix and diagonal

%input = 1:16;

vec = input;
nRows = sqrt(length(vec));
nCols = nRows;

output = zeros(nRows);

r = 1;
c = 1;
currEle = 1;

direction = "right";

while ~(r == nRows && c ==nCols)
    
    switch direction
    
        case "right"
            output(r,c) = vec(currEle);
            c = c+1;
            currEle = currEle+1;
            direction = "digDown";

        case "digDown"
            if c == 1
                output(r,c) = vec(currEle);
                r = r+1;
                currEle = currEle+1;
                direction = "down";
                
            else
                
                output(r,c) = vec(currEle);
                r = r+1;
                c = c-1;
                currEle = currEle+1;
            end
            
        case "down"
            output(r,c) = vec(currEle);
            c = c+1;
            r = r-1;
            currEle = currEle+1;
            direction = "digUp";
    
        case "digUp"
             if r == 1
                direction = "right";
             else
                output(r,c) = vec(currEle);
                r = r-1;
                c = c+1;
                currEle = currEle+1;
             end
    end
    
    if r > nRows
        break
    end
     
end
%length(find(output==0))
%output = output(1:totElements/2);
%end