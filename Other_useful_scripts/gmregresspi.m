function [b,xo,yo,se,pint] = gmregresspi(x,y,xo,alpha)
%GMREGRESSPI Interval Prediction of a Single Value for a Geometric Mean 
%Regression (Reduced Major Axis Regression, RMA).
%   Model II regression should be used when the two variables in the
%   regression equation are random and subject to error, i.e. not 
%   controlled by the researcher. Model I regression using ordinary least 
%   squares underestimates the slope of the linear relationship between the
%   variables when they both contain error. According to Sokal and Rohlf 
%   (1995), the subject of Model II regression is one on which research and
%   controversy are continuing and definitive recommendations are difficult
%   to make.
%
%   GMREGRESSPI is a Model II procedure. It standardize variables before 
%   the slope is computed. Each of the two variables is transformed to have 
%   a mean of zero and a standard deviation of one. The resulting slope is
%   the geometric mean of the linear regression coefficient of Y on X. 
%   Ricker (1973) coined this term and gives an extensive review of Model
%   II regression. Given a comment of the Jolicoeur and Mosimann (1968) and
%   later continued by McArdle (1988). It is also known as Standard(ized) 
%   Major Axis. In short, the OLS slope is divided by the correlation 
%   coefficient.
%
%   Friedman et al. (2013) states that, when there is no scientific reason 
%   to define one variable as dependent on the other in the usual 
%   regression sense, a method such as RMA may be more appropriate. RMA 
%   fits the line to data which minimizes the sum of the areas of the right
%   triangles which have legs parallel to the x-axis, y-axis, and 
%   hypotenuse on the fitted line. Thus given n data points (xi,yi) and the
%   model yi = a + bxi + ei, RMA minimizes,
%
%           i=1_Sum_n(abs((yi - (a + bxi))*(xi - ((yi - a)/b))))
%
%   Complementarily, we broadly recommend you to review the 2009 Smith's 
%   paper 'Use and Misuse of the Reduced Major Axis for Line-Fitting'.
%
%   [B,XO,YO,SE,PINT] = GMREGRESSPI(X,Y,XO,ALPHA) returns the vector B of 
%   regression coefficients in the linear Model II, XO given predictor 
%   value, YO predicted response value, SE standard error and the PINT 
%   interval prediction.
%
%   GMREGRESSPI treats NaNs in X or Y as missing values, and removes them.
%
%   If nargout is empty, it just gives the interval prediction.
%
%   Syntax: function [b,xo,yo,se,pint] = gmregresspi(x,y,xo,alpha)
%
%   Inputs:
%       x - independent variable vector data
%       y - dependent variable vector data
%      xo - given predictor value
%   alpha - significance level (default = 0.05)
%
%   Outputs:
%       b - regression statistics in the linear Model II
%      xo - given predictor value
%      yo - predicted response value
%      se - standard error
%    pint - interval prediction
%
%   Example. From the Box 14.12 (California fish cabezon [Scorpaenichthys 
%   marmoratus]) of Sokal and Rohlf (1995). We are interested to get the
%   95% interval predicted value for a response value of 20. The data are:
%
%   x=[14,17,24,25,27,33,34,37,40,41,42];
%   y=[61,37,65,69,54,93,87,89,100,90,97];
%
%   Calling on Matlab the function: 
%                [b,xo,se,yo,pint] = gmregresspi(x,y,20)
%
%   Answer is:
%
%   b = 12.1938    2.1194
%
%   xo = 20
%
%   yo = 54.5811
%
%   se =  7.8954
%
%   pint = 36.7204   72.4418
%
%   Created by A. Trujillo-Ortiz, R. Hernandez-Walls and D.S. Vega-Valdez 
%             Facultad de Ciencias Marinas
%             Universidad Autonoma de Baja California
%             Apdo. Postal 453
%             Ensenada, Baja California
%             Mexico.
%             atrujo@uabc.edu.mx
%
%   Copyright (C)  April 5, 2014. 
%
%   --We thank Ali Meshgi (Department of Civil & Environmental Engineering,
%     National University of Singapore, for encourage us to write this
%     m-file.--
%
%   To cite this file, this would be an appropriate format:
%   Trujillo-Ortiz, A., R. Hernandez-Walls and D.S. Vega-Valdez. (2014). 
%      gmregresspi:Interval Prediction of a Single Value for a Geometric 
%      Mean Regression (Reduced Major Axis Regression, RMA).  
%      A MATLAB file. [WWW document]. URL http://www.mathworks.com/
%      matlabcentral/fileexchange/26240-gmregresspi
%    
%   References:
%   Barker, F., Soh, Y. E. and Evans, R. J. (1988), Properties of the
%              geometric mean functional relationship. Biometrics, 
%              44:279-281.
%   Friedman, J., Bohonak, A. J. and Levine, R. A. (2013), When are two
%              pieces better than one:fitting and testing OLS and RMA 
%              regressions. Environmetrics, 24:306-316.
%   Jolicoeur, P. and Mosimann, J. E. (1968), Intervalles de confiance pour
%              la pente de l’axe majeur d’une distribution normale 
%              bidimensionnelle. Biométrie-Praximétrie, 9:121-140.
%   McArdle, B. (1988), The structural relationship:regression in biology.
%              Can. Jour. Zool. 66:2329-2339.
%   Ricker, W. E. (1973), Linear regression in fishery research. J. Fish.
%              Res. Board Can., 30:409-434. 
%   Smith, R. J. (2009), Use and Misuse of the Reduced Major Axis for
%              Line-Fitting. Am. Jour Phys. Anthrop. 140:476–486.
%   Sokal, R. R. and Rohlf, F. J. (1995), Biometry. The principles and
%              practice of the statistics in biologicalreserach. 3rd. ed.
%              New-York:W.H.,Freeman. [Sections 14.13 and 15.7] 
%

if  nargin < 3
    error('gmregresspi:TooFewInputs', ...
          'GMREGRESSPI requires at least two input arguments.');
elseif nargin == 3
    alpha = 0.05;
end

x = x(:); y = y(:);

% Check that matrix (X) and rigth hand side (Y) have compatible dimensions
[n,ncolx] = size(x);
if ~isvector(y)
    error('gmregresspi:InvalidData', 'Y must be a vector.');
elseif numel(y) ~= n
    error('gmregresspi:InvalidData', ...
          'The number of rows in Y must equal the number of rows in X.');
end

% Remove missing values, if any
wasnan = (isnan(y) | any(isnan(x),2));
havenans = any(wasnan);
if havenans
   y(wasnan) = [];
   x(wasnan,:) = [];
   n = length(y);
end

R = corrcoef(x,y);
r = R(1,2); %correlation coefficient
s = r/abs(r); %find sign of the correlation coefficient: this former bug 
              %was efficiently corrected thanks to the valuable suggestions
              %given by Holger Goerlitz and Joel E. Cohen. Yes, a negative
              %slope are always negative!
S = cov(x,y);
SCX = S(1,1)*(n-1);
SCY = S(2,2)*(n-1);
SCP = S(1,2)*(n-1);
v = s*sqrt(SCY/SCX); %slope
u = mean(y)-mean(x)*v; %intercept
b = [u v];

%Interval prediction of a single value
%Statistic fundamentals were extracted from Barker et al. (1988) and 
%Friedman et al. (2013)
X = [ones(n,1) x]; %design matrix
ye = X*b';
sde = sqrt(sum(abs((y-ye).*(x-((y-u)/v))))/(n-2)); %residual standard 
                                                   %deviation
se = sde*sqrt(1+(1/n)+(xo-mean(x))^2/SCX); %standard error
t = tinv(1-(alpha/2),n-2);
yo = u + v*xo; %yo-predicted value
li = yo - t*se;
ls = yo + t*se;
pint = [li ls];

return,