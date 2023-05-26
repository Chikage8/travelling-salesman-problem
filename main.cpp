#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>

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
    int currentLocation = -5;
    std::vector<int> locationVisitOrder;
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

//std::vector<Edge> calculateDistances(std::vector<Coordinate>& coordinates) {
//    int numCoordinates = coordinates.size();
//    // create the vector of edges struct
//    std::vector<Edge> edges;
//    // loop over the coordinates to calculate distance and store it in the vector of edges
//    for (int i = 0; i < numCoordinates; ++i) {
//        for (int j = 0; j < numCoordinates; ++j) {
//            if (i != j) {
//                double distance = calculateDistance(coordinates[i], coordinates[j]);                
//                Edge edge;
//                edge.coord1 = coordinates[i];
//                edge.coord2 = coordinates[j];
//                edge.distance = distance;
//                edges.push_back(edge);
//
//            }
//        }
//    }
//}


int main() {
    // set the filename
    std::string filename = "fi10639.tsp";

    // create a vector of Coordinate structs using the readCoordinatesFromFile function which returns vector of Coordinate structs
    std::vector<Coordinate> coordinates = readCoordinatesFromFile(filename);

    // create a vector of Edge structs using the calculateDistances function with the above instance of vector of Coordinate structs as parameter
    //std::vector<Edge> edges = calculateDistances(coordinates);

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

    std::cout << "on first location" << std::endl;
    Tour myFirstTour;
    Tour visitableLocations;
    myFirstTour.currentLocation = 0;
    myFirstTour.locationVisitOrder.push_back(0);
    
    // visitableLocations filled with the ids of coordinates
    for (int i = 0; i < coordinates.size(); i++)
    {
        if (std::find(myFirstTour.locationVisitOrder.begin(), myFirstTour.locationVisitOrder.end(), i) != myFirstTour.locationVisitOrder.end())
        {
            visitableLocations.locationVisitOrder.push_back(i);
        }


        
        /*bool newLocation = true;
        for (int v = 0; v < myFirstTour.locationVisitOrder.size(); v++)
        {
            if (i == v)
            {
                newLocation = false;
                break;
            }
        }
        if (newLocation)
        {

        }*/
    }
    double totalInverseDistanceCurrentEdges = 0;

    for (int i = 0; i < edgeInverseDistanceMatrix[myFirstTour.currentLocation].size(); i++)
    {
        if (i != myFirstTour.currentLocation)
        {
            totalInverseDistanceCurrentEdges += edgeInverseDistanceMatrix[myFirstTour.currentLocation][i];
        }
    }
    
    std::vector<NeighborProbability> probabiltiesOfCurrentEdges;

    for (int i = 0; i < edgeInverseDistanceMatrix[myFirstTour.currentLocation].size(); i++)
    {
        if (i != myFirstTour.currentLocation)
        {
            NeighborProbability neighborProb;
            neighborProb.id = i;
            neighborProb.probability = edgeInverseDistanceMatrix[myFirstTour.currentLocation][i] / totalInverseDistanceCurrentEdges;
            probabiltiesOfCurrentEdges.push_back(neighborProb);
        }
    }

    for (const auto& neighbor : probabiltiesOfCurrentEdges)
    {
        std::cout << neighbor.id << " " << neighbor.probability << std::endl;
    }

    // Print the values of the vector
   /* for (const auto& probability : probabiltiesOfCurrentEdges) {
        std::cout << probability << std::endl;
    }*/


   
    


    // Display the matrix
    /*for (int i = 0; i < numCoordinates; i++) {
        for (int j = 0; j < numCoordinates; j++) {
            std::cout << edgeMatrix[i][0] << " ";
        }
        std::cout << std::endl;
    }*/

    // prints out the vector of Coordinate structs instance "coordinates"
    /*for (const auto& coord : coordinates) {
        std::cout << "ID: " << coord.id << ", X: " << coord.x << ", Y: " << coord.y << std::endl;
    }*/
    

    return 0;
}
