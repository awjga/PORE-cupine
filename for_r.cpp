#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double sd_c(double x_m, double x_s, double x_n,double y_m, double y_s, int y_n)
{
	double al, var, tmp_sd;
            al=x_n+y_n;

            tmp_sd=al*((x_n-1)*(x_s*x_s)+(y_n-1)*(y_s*y_s))+y_n*x_n*(x_m-y_m)*(x_m-y_m);
            var=tmp_sd/(al*(al-1));
            
            return(sqrt(var));
}

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

		sx=sd_c(mx,sx,nx,my,sy,ny);
		mx=(mx*nx+my*ny)/(nx+ny);
		nx=nx+ny;
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
		mx=(mx*nx+my*ny)/(nx+ny);
		nx=nx+ny;
	}

	return(mx);
}	
			
			
			
			
