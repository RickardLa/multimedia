function output = zigZagTest(input)

% Scans upper  left of matrix and diagonal
    

block = input;
[nRows,nCols] = size(input);
totElements = nRows*nCols;

output=zeros(1,totElements);

r = 1;
c = 1;
currEle = 1;

direction = "right";

while ~(r == nRows && c ==nCols)
    
    switch direction
    
        case "right"
            output(currEle) = block(r,c);
            c = c+1;
            currEle = currEle+1;
            direction = "digDown";

        case "digDown"
            if c == 1
                output(currEle) = block(r,c);
                r = r+1;
                currEle = currEle+1;
                direction = "down";
                
            else
                
                output(currEle) = block(r,c);
                r = r+1;
                c = c-1;
                currEle = currEle+1;
            end
            
        case "down"
            output(currEle) = block(r,c);
            c = c+1;
            r = r-1;
            currEle = currEle+1;
            direction = "digUp";
    
        case "digUp"
             if r == 1
                direction = "right";
             else
                output(currEle) = block(r,c);
                r = r-1;
                c = c+1;
                currEle = currEle+1;
             end
    end
    
    if r > nRows
        break
    end
     
end
length(find(output==0))
output = output(1:totElements/2);
end