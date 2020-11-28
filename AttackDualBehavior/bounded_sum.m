function [sum_res] = bounded_sum(vector, from, to)


for i= 1:length(vector)
    
    if (i >=from && i<= to)
        
        vect(i-from+1) = vector(i); 

    end
end

sum_res = sum(vect);

end