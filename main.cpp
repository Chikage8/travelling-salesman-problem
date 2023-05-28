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
    double distanceTraveled = 0;
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

int ChooseStartingCity()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 10638);

    int randomNum = dist(gen);
    return randomNum;
}

int main() {
    // set the filename
    std::string filename = "fi10639.tsp";

    // create a vector of Coordinate structs using the readCoordinatesFromFile function which returns vector of Coordinate structs
    std::vector<Coordinate> coordinates = readCoordinatesFromFile(filename); 

    int numCoordinates = coordinates.size();

    // Create a 2D vector to represent the matrices of inversedistance and edge rating which will be used for choosing neighbors
    std::vector<std::vector<double>> edgeInverseDistanceMatrix(numCoordinates, std::vector<double>(numCoordinates));
    std::vector<std::vector<double>> edgeRatingMatrix(numCoordinates, std::vector<double>(numCoordinates));

    // Fill the inverse distance matrix with distances
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
    
    // Initialize the rating matrix with 1's
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

    double shortestDistance = 2 * std::pow(10, 7);

    for (int run = 0; run < 100; run++)
    {
        // Tour is a struct that contains int(currnetCity), vector of ints(list of visited city ids in order in this case) and a double holding the total distance traveled in the tour
        Tour currentTour;

        // creating set for unvisited cities so that we can remove by value
        std::set<int> visitableCities;

        for (int i = 0; i < coordinates.size(); i++)
        {
            visitableCities.insert(i);
        }

        // choose random starting city
        int startCityIndex = ChooseStartingCity();
        VisitCity(currentTour, visitableCities, startCityIndex);

        while (visitableCities.size() > 0)
        {
            // variable to add all the inverse distances of visitableCities for probability calculation
            double totalInverseDistanceVisitableCities = 0;
            double totalRatingVisitableCities = 0;

            // adding up all the inverse distance for the visitableCities
            for (int city : visitableCities)
            {
                totalInverseDistanceVisitableCities += edgeInverseDistanceMatrix[currentTour.currentCity][city];
                totalRatingVisitableCities += edgeRatingMatrix[city][currentTour.currentCity];
            }

            // vector of NeighborProbability construct to store the id and the probability of cities 
            std::vector<NeighborProbability> probabilityVisitableCities;

            // storing the id and probabilty for each of the visitable cities
            for (int city : visitableCities)
            {
                NeighborProbability neighborProb;
                neighborProb.id = city;
                double neighborRatingProb =  edgeRatingMatrix[currentTour.currentCity][neighborProb.id] / totalRatingVisitableCities;
                double neighborDistanceProb = edgeInverseDistanceMatrix[currentTour.currentCity][neighborProb.id] / totalInverseDistanceVisitableCities;
                neighborProb.probability = neighborDistanceProb * neighborRatingProb;

                probabilityVisitableCities.push_back(neighborProb);
            }
            double totalProbNotNormalized = 0;
            double normalizedTotalProb = 0;
            for (auto& neighbor : probabilityVisitableCities)
            {
                totalProbNotNormalized += neighbor.probability;
            }
            for (auto& neighbor : probabilityVisitableCities)
            {                
                // to get normalized Total Probability
                normalizedTotalProb += neighbor.probability * 1 / totalProbNotNormalized;
                // scale neighbor probability with along Total Probability
                neighbor.probability *= 1 / totalProbNotNormalized;
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
                    currentTour.distanceTraveled += 1 / edgeInverseDistanceMatrix[neighbor.id][currentTour.currentCity];
                    VisitCity(currentTour, visitableCities, neighbor.id);
                    break;
                }
            }
        }
        // Visit back the starting city
        currentTour.distanceTraveled += 1 / edgeInverseDistanceMatrix[startCityIndex][currentTour.currentCity];
        currentTour.currentCity = startCityIndex;
        currentTour.cityVisitOrder.push_back(startCityIndex);

        std::cout << "cityVisitOrder.size() = " << currentTour.cityVisitOrder.size() << std::endl;

        // print current tour distance traveled
        std::cout << "currentTour.distanceTraveled = " << currentTour.distanceTraveled << std::endl;

        // adjusting the edge rating matrix
        for (int i = 0; i < currentTour.cityVisitOrder.size() - 1; i++)
        {
            edgeRatingMatrix[currentTour.cityVisitOrder[i]][currentTour.cityVisitOrder[i + 1]] += std::pow(10, 7) / currentTour.distanceTraveled;
            edgeRatingMatrix[currentTour.cityVisitOrder[i + 1]][currentTour.cityVisitOrder[i]] += std::pow(10, 7) / currentTour.distanceTraveled;
        }

        if (currentTour.distanceTraveled < shortestDistance)
        {
            shortestDistance = currentTour.distanceTraveled;
        }
    }

    return 0;
}







