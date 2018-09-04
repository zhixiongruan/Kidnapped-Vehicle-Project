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

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
  //   x, y, theta and their uncertainties from GPS) and all weights to 1. 
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  
  // Set number of particles
  num_particles = 100;
  
  // This line creates a normal (Gaussian) distribution for x, y, theta
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
   
  // Initializes particles
  for (int i = 0; i < num_particles; ++i) {
    particles[i].id = i;
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
    particles[i].weight = 1.0;
  }
    
  // Initialization completed
  is_initialized = true;


}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/
  
  // creates a normal (Gaussian) distribution
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);
  
  // Add measurements to each particle
  for (int i = 0; i < num_particles; ++i) {
    if (fabs(yaw_rate) < 0.00001) {
	  // avoid divide by zero
      particles[i].x += velocity * delta_t * cos(particles[i].theta);
      particles[i].y += velocity * delta_t * sin(particles[i].theta);
    } else {
      particles[i].x += (velocity/yaw_rate) * (sin(particles[i].theta + (yaw_rate * delta_t)) - sin(particles[i].theta));
      particles[i].y += (velocity/yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate * delta_t)));
      particles[i].theta += yaw_rate * delta_t;
    }

    // Add noise to the particles
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
  
  for (int i=0; i<(int)observations.size(); i++) {
	  double distance = numeric_limits<double>::max(); //the maximum number of double type compiler allowed
	  for (int j=0; j<(int)predicted.size(); j++) {
		double d = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
		if (d<distance) {
			distance = d;
			observations[i].id = predicted[j].id;
		}
	}
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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
  
  /*  Weight update process
	 *  1. Transform the observations into map coodination regarding to the 
	 *  2. Extract the landmark, which are found in sensor range
	 *  3. Associate the observations to the extracted landmarks
	 *  4. Calculate weights
	 */
	
  	for( int i = 0; i < num_particles; i++)
	{
		double x_p = particles[i].x;
		double y_p = particles[i].y;
		double theta_p = particles[i].theta;
		
		/* Transform the observations into map coodination regarding to the*/
  		vector<LandmarkObs> TransformedObs;
		LandmarkObs TransformObs;
		
		for(uint j = 0; j < observations.size(); j++)
		{
			TransformObs.x = cos(theta_p) * observations[j].x - sin(theta_p) * observations[j].y + x_p;
			TransformObs.y = sin(theta_p) * observations[j].x + cos(theta_p) * observations[j].y + y_p;
			TransformObs.id = observations[j].id;
			TransformedObs.push_back(TransformObs);
		}
		
		/* Extract the landmark, which are found in sensor range*/
   		vector<LandmarkObs> ExtractedLandmarks;
		LandmarkObs Landmark;
		for(uint j = 0; j < map_landmarks.landmark_list.size(); j++)
		{
			double Dist = dist(x_p, y_p, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f);
			
			if( Dist < sensor_range)
			{
				Landmark.x = map_landmarks.landmark_list[j].x_f;
				Landmark.y = map_landmarks.landmark_list[j].y_f;
				Landmark.id = map_landmarks.landmark_list[j].id_i;
				
				ExtractedLandmarks.push_back(Landmark);
			}
		}
		
		/* Associate the observations to the extracted landmarks*/
		dataAssociation(ExtractedLandmarks, TransformedObs);
		particles[i].weight = 1.0;
		
		/* Calculate weights */
		for( uint j = 0; j < TransformedObs.size(); j++)
		{
			double x = TransformedObs[j].x;
			double y = TransformedObs[j].y;
			double M_x = 0;
			double M_y = 0;
			/* Looking for the nearest landmark from the extracted landmark*/
			for (uint k = 0; k < ExtractedLandmarks.size(); k++)
			{
				if(TransformedObs[j].id == ExtractedLandmarks[k].id)
				{
					M_x = ExtractedLandmarks[k].x;
					M_y = ExtractedLandmarks[k].y;
				}
			}
			
			particles[i].weight *= exp(-pow(x - M_x,2) / 2 / pow(std_landmark[0], 2) - pow(y - M_y, 2) / 2 / pow(std_landmark[1],2))
									/ 2 / M_PI / std_landmark[0] / std_landmark[1];
		}
		
	}	
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

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
