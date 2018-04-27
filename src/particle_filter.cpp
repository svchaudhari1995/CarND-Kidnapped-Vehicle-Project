/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <cmath>
#include <limits>

#include "particle_filter.h"
# define M_PI 3.14159265358979323846  /* pi */
using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	is_initialized=true;
	default_random_engine gen;
	// cout<<"Init start";
	double std_x, std_y, std_theta; // Standard deviations for x, y, and theta
	num_particles=5;
	// TODO: Set standard deviations for x, y, and theta
	 std_x = std[0];
	 std_y = std[1];
	 std_theta =std[2];
	 

	// This line creates a normal (Gaussian) distribution for x
	normal_distribution<double> dist_x(x, std_x);
	
	// TODO: Create normal distributions for y and theta
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	
	for (int i = 0; i < num_particles; ++i) {
		double sample_x, sample_y, sample_theta;
		
		// TODO: Sample  and from these normal distrubtions like this: 
		//	 sample_x = dist_x(gen);
		//	 where "gen" is the random engine initialized earlier.
		
		 sample_x = dist_x(gen);
		 sample_y = dist_y(gen);
		 sample_theta = dist_theta(gen);	 
		 
		 Particle p;
		 p.id=i;
		 p.x=sample_x;
		 p.y=sample_y;
		 p.weight=1;
		 p.theta=sample_theta;

		 particles.push_back(p);
	}

	// cout<<"Init end";

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	// std::cout<<"hiii";
	default_random_engine gen;
	
	for(int i=0;i<num_particles;i++)
	{
		
	if (fabs(yaw_rate) < 0.00001) {  
    particles[i].x +=  velocity * delta_t * cos(particles[i].theta);
			particles[i].y +=velocity * delta_t * sin(particles[i].theta);
  }
  else{
		particles[i].x+=(velocity/yaw_rate)*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta));
		particles[i].y+=(velocity/yaw_rate)*(-cos(particles[i].theta+yaw_rate*delta_t)+cos(particles[i].theta));
		
  	}
	  
	  particles[i].theta+=yaw_rate*delta_t;
	  normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
		normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
		normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);
		particles[i].x = dist_x(gen);
		 particles[i].y = dist_y(gen);
		particles[i].theta	= dist_theta(gen);	 
	}


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.




}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variatelandmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	
		for(int i=0;i<particles.size();i++)
	{
		vector<double> lsense_x,lsense_y;
		vector<int> lassociations;
		Particle p=particles[i];
		double weight=1.0;
		for(int j=0;j<observations.size();j++)
		{
			double xm=p.x+cos(p.theta)*observations[j].x-sin(p.theta)*observations[j].y;
			double ym=p.y+sin(p.theta)*observations[j].x+cos(p.theta)*observations[j].y;
			double idm=observations[j].id;
			lsense_x.push_back(xm);
			lsense_y.push_back(ym);
			double x_min=0;
			double y_min=0;

			double min=std::numeric_limits<double>::max();
			int minindex;
			for(int k=0;k<map_landmarks.landmark_list.size();k++)
			{
				double d=dist(xm,ym,map_landmarks.landmark_list[k].x_f,map_landmarks.landmark_list[k].y_f);
				
				if(d<min)
				{
					min=d;
					x_min=map_landmarks.landmark_list[k].x_f;
					y_min=map_landmarks.landmark_list[k].y_f;
					minindex=map_landmarks.landmark_list[k].id_i;
				}
	
			}
			double sig_x=std_landmark[0];
			double sig_y=std_landmark[1];
			lassociations.push_back(minindex);
			double gauss_norm= (1/(2 * M_PI * sig_x * sig_y));
			double exponent= pow(xm -x_min,2) /(2 * pow(sig_x,2)) + pow(ym - y_min,2)/(2 * pow(sig_y,2));
			double EulerConstant = std::exp(-exponent);
			weight*= gauss_norm * EulerConstant;
		}
		particles[i].weight=weight;
		particles[i].sense_x=lsense_x;
		particles[i].sense_y=lsense_y;
		particles[i].associations=lassociations;


	}




}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution


	std::random_device rd;     // only used once to initialise (seed) engine
	std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
	double max_alpha=std::numeric_limits<double>::min();
	double total_weight=0.0;
	double local_alpha;
	vector<double> alphas;
	for(int i=0;i<particles.size();i++)
	{

			total_weight+=particles[i].weight;
			
		
	}
	for(int i=0;i<particles.size();i++)
	{
		
		 local_alpha=particles[i].weight/total_weight;
		if(local_alpha>max_alpha)
		{
			max_alpha=local_alpha;
		}
		alphas.push_back(local_alpha);
		
	}

	vector<Particle> new_particles;
	
	std::uniform_int_distribution<int> uni(0,particles.size()-1);
	std::uniform_real_distribution<double> uni_alpha(0.0,2*max_alpha);
	int index=uni(rng);
	double beta=0.0;


	for(int i=0;i<particles.size();i++)
	{
		beta+=uni_alpha(rng);
		while(alphas[index]< beta){
			beta -= alphas[index];
			index++;
			index %= num_particles;
		}
		new_particles.push_back(particles[index]);
	}

	particles = new_particles;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
