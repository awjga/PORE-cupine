#include <Rcpp.h>
using namespace Rcpp;

		
// [[Rcpp::export]]
double sd_combine(NumericVector sd,NumericVector mean, NumericVector freq)
{
	double i , x1, nx, ny, sx, sy, mx, my;
	x1 = mean.size();

	if(x1==1){
		sx=sd[0];
		return(sx);
	}
	mx=mean[0];
	sx=sd[0];
	nx=freq[0];
	for (i =1;i<x1;i++)
	{
		my=mean[i];
		sy=sd[i];
		ny=freq[i];

		sx=sd_combine(mx,sx,nx,my,sy,ny);
		nx=nx+ny;
		mx=(mx*nx+my*ny)/(nx+ny);
	}

	return(sx);
}	

// [[Rcpp::export]]

double mean_combine(NumericVector mean, NumericVector freq)
{
	double i , x1, nx, ny, mx, my;
	x1 = mean.size();

	if(x1==1){
		nx=mean[0];
		return(nx);
	}
	mx=mean[0];
	nx=freq[0];
	for (i =1;i<x1;i++)
	{
		my=mean[i];
		ny=freq[i];

		nx=nx+ny;
		mx=(mx*nx+my*ny)/(nx+ny);
	}

	return(mx);
}	
			
			
			
			