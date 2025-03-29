function zadoffChuSequence = zadoffChu(root, numElement)
% zadoffChu generate root Zadoff-Chu sequence of complex symbol 
% allocated along Frequency Domain
%
%   zadoffChuSequence = zadoffChu(root, numElement) generates the root'th  
%   (physical root sequence index (u)) Zadoff-Chu sequence of length
%   numElement. 
%   
%   Input: root and numElement must be relatively prime    
%   Output: zadoffChuSequence is an numElement-length row
%   vectors of complex symbols.

    if gcd(root, numElement) ~= 1
        warning('GCD(%d, %d) must equal one\n', root, numElement);
        return;
    end
    n = 0:numElement - 1;
    zadoffChuSequence = exp(-1j * ((pi * root * n .* (n + 1)) / numElement));
end