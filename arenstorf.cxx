/**
 * Homework 7
 * Author: Philipp Schulte
 * Date: 14.12.2015
 * 
 * Description:
 * Consider a light body which is moving in the gravitational field of two heavy objects. The motion of the light object will not influence the heavy objects. The heavy bodies with mass ratio 
 * circle in the (x,y) plane with frequency 1 around their common center of gravity, which we assume to be the origin.
 * 
 * The equations of motion of the light body are (see github)
 * 
 * Consider the motion of a light body in the field of earth and moon. The mass ratio is . We restrict ourself to planar motion (z=0) and may use the initial conditions x(0) = 0.994, y(0) = 0,
 * x'(0) = 0 , and y'(0) = -2.00158510637908
 * 
 * This (according to literature) should lead to a periodic motion with period T=17.065216560157. These orbits are known as Arenstorf orbits.
 * 
 * Implement the Dormand-Prince 4/5 RK method with the coefficients (see github)
 * 
 * 
 * to solve the ODEs. The first line of b-coefficients corresponds to the fifth-order result and the second line to the fourth-order result.
 * 
 * Use the norm , where , to implement a step-size control. Limit  to reasonable numbers Tol, i.e. try values in the range of Tol = 1e-3 to Tol = 1e-10. After adjusting the time step continue 
 * the integration using the result from the fourth-order scheme (Check for yourself whether it matters which solution you use to continue the integration).
 * 
 * For Tol = 1e-5 supply a plot for the trajectory (x(t),y(t)) and plot dt versus t to see how the algorithm adjusts the step-size.
 * 
 * 
 */

// Imports and name space description
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>

using namespace std;


// Section to define the needed functions
void func(double* const ki, const double my, const double x, const double x1, const double y, const double y1 ); 
void RungeKytta(const double* f, const double mu, const double dt, double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double* k7);


// main function
int main() {
	// Define variables for calculation
	//First basic constant variables
	const int dim = 4;  
	const double tend = 17.065216560157, epsilon = 1e-5, my = 0.012277471; 
	const double b41=5179.0/57600.0, b42=0.0, b43=7571.0/16695.0, b44=393.0/640.0, b45=-92097.0/339200.0, b46=187.0/2100.0, b47=1.0/40.0;
	const double b51=35.0/384.0, b52=0.0, b53=500.0/1113.0, b54=125.0/192.0, b55=-2187.0/6784.0, b56=11.0/84.0;
	
	
	// the array which contains the calculation results
	
	double R4[dim] = {0.994, 0.0, 0.0, -2.00158510637908};  //Structure of R4 = {x, x', y, y'}
	double R5[dim] = {0.994, 0.0, 0.0, -2.00158510637908};  //Structure of R5 = {x, x', y, y'}
	
	// Helper variables for calculation
	double k1[dim], k2[dim], k3[dim], k4[dim], k5[dim], k6[dim], k7[dim], etemp = 0.0;
	double dt = 1e-3, error, curt = 0;
	int i;
					 
	
	ofstream out("solution");
					 
	
	// Save the values in a file
	out << curt << "\t" << dt << "\t" << R4[0] << "\t" << R4[2] << endl;

	
	while (curt <= tend) {
		error = 0.0;
		// Calculate the K values
		RungeKytta(R4, my, dt, k1, k2, k3, k4, k5, k6, k7);
		
		// Calculate the next step values
		for (i = 0; i < dim; i++) {
			R5[i] = R4[i] + dt*(b51*k1[i]+b52*k2[i]+b53*k3[i]+b54*k4[i]+b55*k5[i]+b56*k6[i]);
			R4[i] += dt*(b41*k1[i]+b42*k2[i]+b43*k3[i]+b44*k4[i]+b45*k5[i]+b46*k6[i]+b47*k7[i]);
		}
		
		// calculate the error with the new calculated values
		for (i = 0; i < dim; i++) {
			etemp = abs(R4[i] - R5[i]);
			if (etemp > error) error = etemp;
		}
		
		// Make stepsize controll
		dt *= pow((epsilon/error), 0.2);
		curt += dt;
		out << curt << "\t" << dt << "\t" << R4[0] << "\t" << R4[2] << endl;
	}
	
	out.close();
	
	return 0;
}

// Section to implement the functions

/*
 * 
 * 
 */
void func(double* const ki, const double my, const double x, const double x1, const double y, const double y1 ) {
	
	// First calculate r ans s for this step:
	double r = sqrt(((x + my) * (x + my)) + (y * y));	
	double s = sqrt(((x - 1.0 + my) * (x - 1.0 + my)) + (y * y));	
	
	ki[0] = x1;
	ki[1] = x + 2.0 * y1 - (((1.0 - my) * (x + my)) / (r * r * r)) - ((my * (x - 1.0 + my)) / (s * s * s));
	ki[2] = y1;
	ki[3] = y - 2.0 * x1 - ((y * (1.0 - my)) / (r * r * r)) - ((my * y) / (s * s * s));
}

void RungeKytta(const double* f, const double my, const double dt, double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double* k7) {
	
	// Save the values if my R4 array in variables which are better to identicate.
	double x=f[0], x1=f[1], y=f[2], y1=f[3];
	
	// Set the coefficitents for the calcualtion
	const double a21=1.0/5.0;
	const double a31=3.0/40.0, a32=9.0/40.0;
	const double a41=44.0/45.0, a42=-56.0/15.0, a43=32.0/9.0;
	const double a51=19372.0/6561.0, a52=-25360.0/2187.0, a53=64448.0/6561.0, a54=-212.0/729.0;
	const double a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0, a64=49.0/176.0, a65=-5103.0/18656.0;
	const double a71=35.0/384.0, a72=0.0, a73=500.0/1113.0, a74=125.0/192.0, a75=-2187.0/6784.0, a76=11.0/84.0;
	
	// Do calculation
	func(k1, my, x, x1, y, y1);
	func(k2, my, x+dt*(a21*k1[0]), x1+dt*(a21*k1[1]), y+dt*(a21*k1[2]), y1+dt*(a21*k1[3]));
	func(k3, my, x+dt*(a31*k1[0]+a32*k2[0]), x1+dt*(a31*k1[1]+a32*k2[1]), y+dt*(a31*k1[2]+a32*k2[2]), y1+dt*(a31*k1[3]+a32*k2[3]));
	func(k4, my, x+dt*(a41*k1[0]+a42*k2[0]+a43*k3[0]), x1+dt*(a41*k1[1]+a42*k2[1]+a43*k3[1]), y+dt*(a41*k1[2]+a42*k2[2]+a43*k3[2]), y1+dt*(a41*k1[3]+a42*k2[3]+a43*k3[3]));
	func(k5, my, x+dt*(a51*k1[0]+a52*k2[0]+a53*k3[0]+a54*k4[0]), x1+dt*(a51*k1[1]+a52*k2[1]+a53*k3[1]+a54*k4[1]), y+dt*(a51*k1[2]+a52*k2[2]+a53*k3[2]+a54*k4[2]), y1+dt*(a51*k1[3]+a52*k2[3]+a53*k3[3]+a54*k4[3]));
	func(k6, my, x+dt*(a61*k1[0]+a62*k2[0]+a63*k3[0]+a64*k4[0]+a65*k5[0]), x1+dt*(a61*k1[2]+a62*k2[1]+a63*k3[1]+a64*k4[1]+a65*k5[1]), y+dt*(a61*k1[2]+a62*k2[2]+a63*k3[2]+a64*k4[2]+a65*k5[2]), y1+dt*(a61*k1[3]+a62*k2[3]+a63*k3[3]+a64*k4[3]+a65*k5[3]));
	func(k7, my, x+dt*(a71*k1[0]+a72*k2[0]+a73*k3[0]+a74*k4[0]+a75*k5[0]+a76*k6[0]), x1+dt*(a71*k1[1]+a72*k2[1]+a73*k3[1]+a74*k4[1]+a75*k5[1]+a76*k6[1]), y+dt*(a71*k1[2]+a72*k2[2]+a73*k3[2]+a74*k4[2]+a75*k5[2]+a76*k6[2]), y1+dt*(a71*k1[3]+a72*k2[3]+a73*k3[3]+a74*k4[3]+a75*k5[3]+a76*k6[3]));
	
}

