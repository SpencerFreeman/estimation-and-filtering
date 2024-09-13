%  centrall.m Matlab commands.  MAE 676 homework to test the central limit theorem
%
%
%  These Matlab commands produce 100 different 100,000-point sample realizations
%  of a sinusoidally-distributed random variable, y.  The column vector
%  y(:,i) contains 100,000 sample realizations of 
%  the sinusoidal distribution.  The column vector theta(:,i)
%  contains the randomized phases (randomized using a flat distribution) that
%  are used to generate the y(:,i) samples.  The column vector
% 
%           z = (1/sqrt(100))(y(:,1) + y(:,2) + ... + y(:,100))
%  
%  contains 100,000 points of a random variable that is the average of
%  100 independent sinusoidally distributed random variables.
%  
%  Execute these lines and plot the results in order to see what a sinusoidal
%  distribution looks like and to see whether the distribution of z approaches
%  a Gaussian, as predicted by the central limit theorem.
%  Plot the histogram of y(:,1), as shown below, and plot the histogram
%  for some other column of y.  Also, plot the histogram of z.  Hand in all
%  3 histograms.  Always add axis labels and titles to you plots.
%
%  Note: the hist function essentially plots the probability density function,
%  except that the vertical scale's units are in number of counts rather than in
%  units of probability per unit of z.  See if you can transform into the latter
%  probability density units.
%
%
» theta = 2*pi*rand(100000,100);
» y = sin(theta);
» hist(y(:,1),20)
» z = sum(y,2)*(1/sqrt(100));
» hist(z,20)