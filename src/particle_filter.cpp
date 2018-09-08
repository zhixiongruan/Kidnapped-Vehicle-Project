/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

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
	for (int i = 0; i < num_particles; i++) {
		Particle particle;
		particle.id = i;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = 1.0;
		particles.push_back(particle);
		weights.push_back(1.0);
	}

	// Initialization completed
	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[],
		double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// Creates a normal (Gaussian) distribution
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);

	// Add measurements to each particle
	for (int i = 0; i < num_particles; i++) {
		if (fabs(yaw_rate) < 0.00001) {
			// Avoid divide by zero
			particles[i].x += velocity * delta_t * cos(particles[i].theta);
			particles[i].y += velocity * delta_t * sin(particles[i].theta);
		} else {
			particles[i].x += (velocity / yaw_rate)
					* (sin(particles[i].theta + (yaw_rate * delta_t))
							- sin(particles[i].theta));
			particles[i].y += (velocity / yaw_rate)
					* (cos(particles[i].theta)
							- cos(particles[i].theta + (yaw_rate * delta_t)));
			particles[i].theta += yaw_rate * delta_t;
		}

		// Add noise to the particles
		particles[i].x += dist_x(gen);
		particles[i].y += dist_y(gen);
		particles[i].theta += dist_theta(gen);

	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted,
		std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	for (int i = 0; i < (int) observations.size(); i++) {
		double distance = numeric_limits<double>::max(); //the maximum number of double type compiler allowed
		for (int j = 0; j < (int) predicted.size(); j++) {
			double d = dist(observations[i].x, observations[i].y,
					predicted[j].x, predicted[j].y);
			if (d < distance) {
				distance = d;
				observations[i].id = predicted[j].id;
			}
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations,
		const Map &map_landmarks) {
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

	double sig_x_pow = pow(std_landmark[0], 2);
	double sig_y_pow = pow(std_landmark[1], 2);
	double gauss_norm = 1.0 / (2.0 * M_PI * std_landmark[0] * std_landmark[1]);
	double weight_sum = 0.0;

	for (int i = 0; i < num_particles; i++) {
		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta;

		// Transformed into map coordinates
		vector<LandmarkObs> map_observations;
		for (int j = 0; j < (int) observations.size(); j++) {
			LandmarkObs lo;
			//Transform to map x coordinate
			lo.x = x + cos(theta) * observations[j].x
					- sin(theta) * observations[j].y;
			//Transform to map y coordinate
			lo.y = y + sin(theta) * observations[j].x
					+ cos(theta) * observations[j].y;
			lo.id = observations[j].id;
			map_observations.push_back(lo);
		}

		//Search for landmarks in particle range
		vector<LandmarkObs> lo_predicted;
		for (int j = 0; j < (int) map_landmarks.landmark_list.size(); j++) {
			double d = dist(x, y, map_landmarks.landmark_list[j].x_f,
					map_landmarks.landmark_list[j].y_f);
			if (d < sensor_range) {
				LandmarkObs lo;
				lo.x = map_landmarks.landmark_list[j].x_f;
				lo.y = map_landmarks.landmark_list[j].y_f;
				lo.id = map_landmarks.landmark_list[j].id_i;
				lo_predicted.push_back(lo);
			}
		}

		// Associate the observations to the landmarks
		dataAssociation(lo_predicted, map_observations);

		particles[i].weight = 1.0;
		// Calculate weights
		for (int j = 0; j < (int)map_observations.size(); j++) {
			for (int k = 0 ; k < (int)lo_predicted.size(); k++){
				if (map_observations[j].id == lo_predicted[k].id) {
					//calculate exponent
					double exponent = (pow((map_observations[j].x - lo_predicted[k].x), 2)) / (2.0 * sig_x_pow)
										+ (pow((map_observations[j].y - lo_predicted[k].y), 2)) / (2.0 * sig_y_pow);
					//calculate weight using normalization terms and exponent
					particles[i].weight *= gauss_norm * exp(-1.0*exponent);
				}
			}
		}
		weight_sum += particles[i].weight;
	}

	// Normalize weights
	for (int i = 0; i < (int)particles.size(); i++) {
		particles[i].weight /= weight_sum;
		weights[i] = particles[i].weight;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	// use discrete_distribution to do the resampling
	discrete_distribution<int> dist(weights.begin(), weights.end());
	vector<Particle> new_particles;
	for (int i = 0; i < num_particles; i++) {
		new_particles.push_back(particles[dist(gen)]);
	}
	particles = new_particles;

}

Particle ParticleFilter::SetAssociations(Particle& particle,
		const std::vector<int>& associations,
		const std::vector<double>& sense_x,
		const std::vector<double>& sense_y) {
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations = associations;
	particle.sense_x = sense_x;
	particle.sense_y = sense_y;

	return particle;
}

string ParticleFilter::getAssociations(Particle best) {
	vector<int> v = best.associations;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1);  // get rid of the trailing space
	return s;
}
string ParticleFilter::getSenseX(Particle best) {
	vector<double> v = best.sense_x;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1);  // get rid of the trailing space
	return s;
}
string ParticleFilter::getSenseY(Particle best) {
	vector<double> v = best.sense_y;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1);  // get rid of the trailing space
	return s;
}
