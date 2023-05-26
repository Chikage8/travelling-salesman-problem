#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include "main.h"
#include <set>

// defining a struct for storing coordinates
struct Coordinate {
    int id;
    double x;
    double y;
};

// defining a struct for storing edges
struct Edge {
    Coordinate coord1;
    Coordinate coord2;
    double distance;
};

struct Tour {    
    int currentCity = -5;
    std::vector<int> cityVisitOrder;
};

struct NeighborProbability {
    int id;
    double probability;
};

std::vector<Coordinate> readCoordinatesFromFile(const std::string& filename) {

    // read file
    std::ifstream file(filename);

    // create the vector of coordinates struct
    std::vector<Coordinate> coordinates;

    if (file.is_open()) {
        std::string line;
        bool startReading = false;

        while (std::getline(file, line)) {
            // when the NODE_COORD_SECTION line is met, set startReading as true and continue with the next line where coordinates actually begin
            if (line.find("NODE_COORD_SECTION") != std::string::npos) {
                startReading = true;
                continue;
            }

            // start reading the coords when the conditions are met
            if (startReading && !line.empty()) {
                // istringstream is used to get to string input and create a Coordinate struct instance out of it to store it later on the coordinates vector of Coordinate structs
                std::istringstream iss(line);
                Coordinate coord;
                if (iss >> coord.id >> coord.x >> coord.y) {
                    coord.id -= 1;
                    coordinates.push_back(coord);
                }
            }
        }

        file.close();
    }
    else {
        std::cerr << "Unable to open the file: " << filename << std::endl;
    }

    return coordinates;
}

// just returns the distance between two coordinate struct instances
double calculateDistance(const Coordinate& coord1, const Coordinate& coord2) {
    double xDiff = coord2.x - coord1.x;
    double yDiff = coord2.y - coord1.y;
    double distance = std::sqrt(xDiff * xDiff + yDiff * yDiff);
    return distance;
}

std::vector<Edge> calculateDistances(std::vector<Coordinate>& coordinates) {
    int numCoordinates = coordinates.size();
    // create the vector of edges struct
    std::vector<Edge> edges;
    // loop over the coordinates to calculate distance and store it in the vector of edges
    for (int i = 0; i < numCoordinates; ++i) {
        for (int j = 0; j < numCoordinates; ++j) {
            if (i != j) {
                double distance = calculateDistance(coordinates[i], coordinates[j]);                
                Edge edge;
                edge.coord1 = coordinates[i];
                edge.coord2 = coordinates[j];
                edge.distance = distance;
                edges.push_back(edge);

            }
        }
    }
    return edges;
}

double generateRandomProb()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // Generate a random number between 0 and 1
    double randomNumber = dis(gen);

    return randomNumber;
}

bool foundIn(Tour& currentTour, const int& i)
{
    return std::find(currentTour.cityVisitOrder.begin(), currentTour.cityVisitOrder.end(), i) != currentTour.cityVisitOrder.end();
}

void VisitCity(Tour& currentTour, std::set<int>& visitableCities, int i)
{
    currentTour.currentCity = i;
    currentTour.cityVisitOrder.push_back(i);

    visitableCities.erase(i);
}


int main() {
    // set the filename
    std::string filename = "fi10639.tsp";

    // create a vector of Coordinate structs using the readCoordinatesFromFile function which returns vector of Coordinate structs
    std::vector<Coordinate> coordinates = readCoordinatesFromFile(filename); 

    int numCoordinates = coordinates.size();

    // Create a 2D vector to represent the matrix
    std::vector<std::vector<double>> edgeInverseDistanceMatrix(numCoordinates, std::vector<double>(numCoordinates));
    std::vector<std::vector<double>> edgeRatingMatrix(numCoordinates, std::vector<double>(numCoordinates));

    // Populate the matrix with distances
    for (int i = 0; i < numCoordinates; i++) {
        for (int j = 0; j < numCoordinates; j++) {
            if (i == j) {
                edgeInverseDistanceMatrix[i][j] = 0.0;  // Distance from a coordinate to itself is 0
            }
            else {
                edgeInverseDistanceMatrix[i][j] = 1 / calculateDistance(coordinates[i], coordinates[j]);
            }
        }
    }
    
    for (int i = 0; i < numCoordinates; i++) {
        for (int j = 0; j < numCoordinates; j++) {
            if (i == j) {
                edgeRatingMatrix[i][j] = 0.0;  // Distance from a coordinate to itself is 0
            }
            else {
                edgeRatingMatrix[i][j] = 1;
            }
        }
    }    

    // my first tour to get the algorithm going 
    // Tour is a struct that contains int(currnetCity) and vector of ints(list of visited city ids in order in this case)
    Tour currentTour;

    // creating set for unvisited cities so that we can remove by value
    std::set<int> visitableCities;

    for (int i = 0; i < coordinates.size(); i++)
    {
        visitableCities.insert(i);        
    }

    VisitCity(currentTour, visitableCities, 0);
    
    while (visitableCities.size() > 0)
    {
        // variable to add all the inverse distances of visitableCities for probability calculation
        double totalInverseDistanceVisitableCities = 0;

        // adding up all the inverse distance for the visitableCities
        for (int city : visitableCities)
        {
            totalInverseDistanceVisitableCities += edgeInverseDistanceMatrix[currentTour.currentCity][city];
        }     
    
        // vector of NeighborProbability construct to store the id and the probability of cities 
        std::vector<NeighborProbability> probabilityVisitableCities;

        // storing the id and probabilty for each of the visitable cities
        for (int city : visitableCities)
        {
            NeighborProbability neighborProb;
            neighborProb.id = city;
            neighborProb.probability = edgeInverseDistanceMatrix[currentTour.currentCity][city] / totalInverseDistanceVisitableCities;

            probabilityVisitableCities.push_back(neighborProb);
        }       

        // Generating a random number between 0 and 1
        double roll = generateRandomProb();

        // declaring a variable to hold the sum of probabilities
        double cumulativeProbability = 0;

        // Choose City based on the probabilty of each neighbor
        // Adding up the probabilty of each visitable city until we exceed the roll
        // When we have exceeded the roll the city that caused cumulativeProbability to exceed the roll is choosed as the next city to visit
        for (const NeighborProbability& neighbor : probabilityVisitableCities)
        {
            cumulativeProbability += neighbor.probability;
            if (cumulativeProbability > roll)
            {
                VisitCity(currentTour, visitableCities, neighbor.id);
                break;
            }            
        }    
        
    }
    for (int i : currentTour.cityVisitOrder)
    {
        std::cout << i << std::endl;
    }

    return 0;
}





