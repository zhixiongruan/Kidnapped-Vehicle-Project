#include <iostream>
#include <algorithm>
#include <vector>

#include "helpers.h"
using namespace std;

std::vector<float> initialize_priors(int map_size, std::vector<float> landmark_positions,
                                     float position_stdev);

float motion_model(float pseudo_position, float movement, std::vector<float> priors,
                   int map_size, int control_stdev);

int main() {
    
    //set standard deviation of control:
    float control_stdev = 1.0f;
    
    //set standard deviation of position:
    float position_stdev = 1.0f;

    //meters vehicle moves per time step
    float movement_per_timestep = 1.0f;

    //number of x positions on map
    int map_size = 25;

    //initialize landmarks
    std::vector<float> landmark_positions {5, 10, 20};
    
    // initialize priors
    std::vector<float> priors = initialize_priors(map_size, landmark_positions,
                                                  position_stdev);
    
    //step through each pseudo position x (i)    
    for (unsigned int i = 0; i < map_size; ++i) {
        float pseudo_position = float(i);

        //get the motion model probability for each x position
        float motion_prob = motion_model(pseudo_position, movement_per_timestep,
                            priors, map_size, control_stdev);
        
        //print to stdout
        std::cout << pseudo_position << "\t" << motion_prob << endl;
    }    

    return 0;
};

//motion model: calculates prob of being at an estimated position at time t
float motion_model(float pseudo_position, float movement, std::vector<float> priors,
                   int map_size, int control_stdev) {

    //initialize probability
    float position_prob = 0.0f;

    //loop over state space for all possible positions x (convolution):
    for (unsigned int j=0; j< map_size; ++j) {
        float next_pseudo_position = float(j);
        //distance from i to j
        float distance_ij = pseudo_position-next_pseudo_position;

        //transition probabilities:
        float transition_prob = Helpers::normpdf(distance_ij, movement, 
                            control_stdev);
        //estimate probability for the motion model, this is our prior
        position_prob += transition_prob*priors[j];
    }
    return position_prob;
}

//initialize priors assumimg vehicle at landmark +/- 1.0 meters position stdev
std::vector<float> initialize_priors(int map_size, std::vector<float> landmark_positions,
                                     float position_stdev) {
//initialize priors assumimg vehicle at landmark +/- 1.0 meters position stdev

    //set all priors to 0.0
    std::vector<float> priors(map_size, 0.0);

    //set each landmark positon +/1 to 1.0/9.0 (9 possible postions)
    float normalization_term = landmark_positions.size() * (position_stdev * 2 + 1);
    for (unsigned int i = 0; i < landmark_positions.size(); i++){
        int landmark_center = landmark_positions[i];
        priors[landmark_center] = 1.0f/normalization_term;
        priors[landmark_center - 1] = 1.0f/normalization_term;
        priors[landmark_center + 1] = 1.0f/normalization_term;

    }
    return priors;
}


#include <iostream>
#include <algorithm>
#include <vector>

#include "helpers.h"
using namespace std;

//set standard deviation of control:
float control_stdev = 1.0f;

//meters vehicle moves per time step
float movement_per_timestep = 1.0f;

//number of x positions on map
int map_size = 25;

//define landmarks
std::vector<float> landmark_positions {5, 10, 12, 20};



std::vector<float> pseudo_range_estimator(std::vector<float> landmark_positions, float pseudo_position);


int main() {        

    //step through each pseudo position x (i)
    for (unsigned int i = 0; i < map_size; ++i) {
        float pseudo_position = float(i);
        //get pseudo ranges
        std::vector<float> pseudo_ranges = pseudo_range_estimator(landmark_positions, pseudo_position);

        //print to stdout
        if (pseudo_ranges.size() >0) {
            for (unsigned int s = 0; s < pseudo_ranges.size(); ++s) {
                std::cout << "x: " << i << "\t" << pseudo_ranges[s] << endl;
            }
            std::cout << "-----------------------" << endl;
        }   
    } 

    return 0;
};

std::vector<float> pseudo_range_estimator(std::vector<float> landmark_positions, float pseudo_position) {
    
    //define pseudo observation vector:
    std::vector<float> pseudo_ranges;
            
    //loop over number of landmarks and estimate pseudo ranges:
        for (unsigned int l=0; l< landmark_positions.size(); ++l) {

            //estimate pseudo range for each single landmark 
            //and the current state position pose_i:
            float range_l = landmark_positions[l] - pseudo_position;
            
            //check if distances are positive: 
            if (range_l > 0.0f) {
                pseudo_ranges.push_back(range_l);
            }
        }

    //sort pseudo range vector:
    sort(pseudo_ranges.begin(), pseudo_ranges.end());
    
    return pseudo_ranges;
}


#include <iostream>
#include <algorithm>
#include <vector>

#include "helpers.h"
using namespace std;

//function to get pseudo ranges
std::vector<float> pseudo_range_estimator(std::vector<float> landmark_positions, 
                                          float pseudo_position);

//observation model: calculates likelihood prob term based on landmark proximity
float observation_model(std::vector<float> landmark_positions, std::vector<float> observations, 
                        std::vector<float> pseudo_ranges, float distance_max, 
                        float observation_stdev);


int main() {  

    //set observation standard deviation:
    float observation_stdev = 1.0f;

    //number of x positions on map
    int map_size = 25;

    //set distance max
    float distance_max = map_size;

    //define landmarks
    std::vector<float> landmark_positions {5, 10, 12, 20};

    //define observations
    std::vector<float> observations {5.5, 13, 15};

    //step through each pseudo position x (i)
    for (unsigned int i = 0; i < map_size; ++i) {
        float pseudo_position = float(i);

        //get pseudo ranges
        std::vector<float> pseudo_ranges = pseudo_range_estimator(landmark_positions, 
                                                                  pseudo_position);

        //get observation probability
        float observation_prob = observation_model(landmark_positions, observations, 
                                                   pseudo_ranges, distance_max, 
                                                   observation_stdev);

        //print to stdout
        std::cout << observation_prob << endl; 
    }      

    return 0;
};

//observation model: calculates likelihood prob term based on landmark proximity
float observation_model(std::vector<float> landmark_positions, std::vector<float> observations, 
                        std::vector<float> pseudo_ranges, float distance_max,
                        float observation_stdev) {

    //initialize observation probability:
    float distance_prob = 1.0f;

    //run over current observation vector:
    for (unsigned int z=0; z< observations.size(); ++z) {

        //define min distance:
        float pseudo_range_min;
        
        //check, if distance vector exists:
        if(pseudo_ranges.size() > 0) {
            //set min distance:
            pseudo_range_min = pseudo_ranges[0];
            //remove this entry from pseudo_ranges-vector:
            pseudo_ranges.erase(pseudo_ranges.begin());

        }    

    //no or negative distances: set min distance to a large number:
    else {

        pseudo_range_min = std::numeric_limits<const float>::infinity();

    }

        //estimate the probabiity for observation model, this is our likelihood: 
        distance_prob *= Helpers::normpdf(observations[z], pseudo_range_min,
                                          observation_stdev);
       
    }
    return distance_prob;
}

std::vector<float> pseudo_range_estimator(std::vector<float> landmark_positions,
                                          float pseudo_position) {
    
    //define pseudo observation vector:
    std::vector<float> pseudo_ranges;
            
    //loop over number of landmarks and estimate pseudo ranges:
        for (unsigned int l=0; l< landmark_positions.size(); ++l) {

            //estimate pseudo range for each single landmark 
            //and the current state position pose_i:
            float range_l = landmark_positions[l] - pseudo_position;
            
            //check if distances are positive: 
            if (range_l > 0.0f) {
                pseudo_ranges.push_back(range_l);
            }
        }

    //sort pseudo range vector:
    sort(pseudo_ranges.begin(), pseudo_ranges.end());
    
    return pseudo_ranges;
}