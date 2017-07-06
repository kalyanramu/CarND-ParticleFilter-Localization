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
#include "particle_filter.h"

#define N_PARTICLES 100

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	std::normal_distribution<double> x_dist(x,std[0]);
	std::normal_distribution<double> y_dist(y,std[0]);
	std::normal_distribution<double> theta_dist(theta,std[0]);
	std::default_random_engine gen;

	num_particles = N_PARTICLES;
	cout << "Initializing Particles" << endl;
	//cout << "Number of particles: "<< particles.size() << endl;
	particles.resize(num_particles);
	weights.resize(num_particles);

	for (int i=0; i < num_particles; i++)
	{
		particles[i].x = x + x_dist(gen) ;
		particles[i].y = y + y_dist(gen);
		particles[i].theta = theta + theta_dist(gen);
		//cout << "particle_" << i << ".x= " << particles[i].x << endl;
	}

	//End of Initialization
	cout << "End of Initialization" << endl;
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	//Velocity and Yaw rate are control inputs measured using motion sensors
	float v = velocity;
	float theta_dot = yaw_rate;
	float theta0,theta1; //state vector

	//Gaussian Noise
	std::normal_distribution<double> dist_x(0,std_pos[0]); 
	std::normal_distribution<double> dist_y(0,std_pos[1]); 
	std::normal_distribution<double> dist_theta(0,std_pos[2]); 
	std::default_random_engine gen;

	//Estimate new state
	for (int i=0; i < num_particles; i++)
	{
		if(yaw_rate !=0)
		{
			theta0 = particles[i].theta;
			theta1 = theta0 + theta_dot * delta_t ;
			particles[i].x = particles[i].x + v/(theta_dot) * (sin(theta1) - sin(theta0));
			particles[i].y = particles[i].y + v/(theta_dot) * (cos(theta0) - cos(theta1));
			particles[i].theta = theta1;
		}
		else
		{
			particles[i].x = particles[i].x + v *delta_t;
			particles[i].y = particles[i].y + v *delta_t;
		}

		//Add gaussian noise to particles
		particles[i].x += dist_x(gen);
		particles[i].y += dist_y(gen);
		particles[i].theta += dist_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
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

	// For each particle :
	// Convert observation/measurements from particle co-ordinate system to map/global co-ordinate system
	// Find associated landmark for each observation/measurement
	// Get the probability of particle taking such measurement using sensor standard deviation & range

	// Learn How to use: sensor range, landmark obsd std
	//sensor range is used to filter out particles beyond sensor's range w.r.t landmarks
	//landmark_std is used for weight calulation of particle as std
	//Step 1: convert sensor measurements to global co-ordinate system
	LandmarkObs global_obs;
	double final_weight,normalized_weight, exp_term;
	double x_diff, y_diff,currobs_weight;
	double x_std_lmx_sqr= std_landmark[0],y_std_lmx_sqr=std_landmark[1];
	double prob_denom = 2*M_PI*std_landmark[0]*std_landmark[1];
	Particle p;
	landmark nearest_landmark;

	for(int i =0; i < num_particles;i++)
	{
		final_weight = 1.0;
		p = particles[i];
		for(LandmarkObs obs: observations)
		{
			//Step 1: Convert observations to global co-ordinate system with reference to particle
			//Fixed bug: sin & cos projections should be done for obs & not for particle positions
			global_obs.x = p.x + (obs.x * cos(p.theta) - obs.y *sin(p.theta));
			global_obs.y = p.y + (obs.x * sin(p.theta) + obs.y *cos(p.theta));

			//Step 2: Find the closest landmark in the map for the particle
			//sensor range is used to filter landmarks that are beyond sensor range from the particle
			//This filtering speeds up computation
			auto nearest_landmark= closestLandMarkPos(global_obs.x, global_obs.y, map_landmarks);

			//Step 3: Compute probablibilty of nearest landmark using multi-variate gaussian
			x_diff = (global_obs.x - nearest_landmark.x);
			y_diff = (global_obs.y - nearest_landmark.y);
			exp_term = (x_diff*x_diff/x_std_lmx_sqr) + (y_diff*y_diff/y_std_lmx_sqr);		
			currobs_weight = exp(-0.5*exp_term);

			//Ste4 : Update the final weight for the particle
			final_weight *= currobs_weight/prob_denom;
		}
		weights[i] = final_weight;
	}

	//Step 5: Normalize weights for all the particles so that cummulative prob =1
	double weights_sum =0;
	for (auto weight: weights)
		weights_sum += weight;

	for (int i =0; i < weights.size();i++)
	{
		
		normalized_weight = weights[i]/ weights_sum;
		weights[i] = normalized_weight;
		particles[i].weight = normalized_weight;
	}

}

landmark ParticleFilter::closestLandMarkPos(double x, double y, Map map_landmarks)
{
	landmark nearest_landmark; int idx;
	int num_landmarks = map_landmarks.landmark_list.size();
	double dist,xm,ym,xdiff,ydiff;
	double nearest_dist = -1.0;

	for (int i=0; i < num_landmarks; i++)
	{
		xm = map_landmarks.landmark_list[i].x_f;
		ym = map_landmarks.landmark_list[i].y_f;
		idx = map_landmarks.landmark_list[i].id_i;
		xdiff = xm - x;
		ydiff = ym - y;
		dist = xdiff*xdiff + ydiff*ydiff;
		//cout << idx << " Dist:" << dist << endl;
		if (dist < nearest_dist || i == 0)
		{
			nearest_landmark.x = xm;
			nearest_landmark.y = ym;
			nearest_landmark.id = idx;
			nearest_dist = dist;
		}

	}
	//cout << "Nearest Landmark: " << nearest_landmark.id << endl;
	return nearest_landmark;
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	//Step 1: Get new particles
    std::discrete_distribution<> weights_cdf(weights.begin(),weights.end());
    std::default_random_engine gen;
	int random_idx;
	std::vector<Particle> new_particles;
	for (int i=0; i < num_particles; i++)
	{
		random_idx = weights_cdf(gen);
		new_particles.push_back(particles[random_idx]);
	}

	//Step 2: Replace old particles with new particles
	particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
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
