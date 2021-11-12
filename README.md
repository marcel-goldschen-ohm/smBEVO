# smBEVO
Computer vision approach to intuitive baseline drift estimation for single-molecule data series.

[![View smBEVO on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/101904-smbevo)

---

### Script

Download `smBEVO.m` and put it in yout MATLAB path. Now you can call the smBEVO function in your scripts.

```matlab
% data ==> Series column data in [x y] or [y] format.
% sigmaX ==> Estimate as longest period in x where baseline is assured to not vary much.
% sigmaY ==> Estimate as ~1/3 of the smallest level separation to be detected in y.
% There are also additional optional arguments, 
%     but setting sigmaX and sigmaY as appropriate for the data is critical.
% results is a struct containing the baseline estimation (i.e. bottom level)
results = smBEVO(data, sigmaX, sigmaY);
```

### MATLAB APP

Install the MATLAB app via the `smBEVO.mlappinstall` file. This provides a graphical interface for loading and baselineing data series.



