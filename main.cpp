#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include "main.h"
#include <set>
#include <queue>
#include <stack>

// defining a struct for storing coordinates
struct Coordinate 
{
    int id;
    double x;
    double y;
};

// defining a struct for storing edges
struct Edge 
{
    Coordinate coord1;
    Coordinate coord2;
    double distance;
};

struct Tour 
{    
    int currentCity = -5;
    double distanceTraveled = 0;
    std::vector<int> cityVisitOrder;
};

struct NeighborProbability 
{
    int id;
    double probability;
};

Tour bestTour;

std::vector<Coordinate> readCoordinatesFromFile(const std::string& filename) 
{

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

// find root in a union find data structure
int findRoot(std::vector<int>& parent, int node)
{
    if (parent[node] != node)
        parent[node] = findRoot(parent, parent[node]);
    return parent[node];
}

// performs union operation
void unionNodes(std::vector<int>& parent, std::vector<int>& rank, int node1, int node2)
{
    int root1 = findRoot(parent, node1);
    int root2 = findRoot(parent, node2);

    if (rank[root1] < rank[root2])
        parent[root1] = root2;
    else if (rank[root1] > rank[root2])
        parent[root2] = root1;
    else
    {
        parent[root2] = root1;
        rank[root1]++;
    }
}

// returns all the possible edges
std::vector<Edge> calculateDistances(std::vector<Coordinate>& coordinates) {
    int numCoordinates = coordinates.size();
    // create the vector of edges struct
    std::vector<Edge> edges;
    // loop over the coordinates to calculate distance and store it in the vector of edges
    for (int i = 0; i < numCoordinates; ++i) {
        for (int j = i+1; j < numCoordinates; ++j) {            
                double distance = calculateDistance(coordinates[i], coordinates[j]);                
                Edge edge;
                edge.coord1 = coordinates[i];
                edge.coord2 = coordinates[j];
                edge.distance = distance;
                edges.push_back(edge);
        }
    }
    return edges;
}

// used to generate random number between 0 and 1
double generateRandomProb()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // Generate a random number between 0 and 1
    double randomNumber = dis(gen);

    return randomNumber;
}

// to find whether the city(i) is present in the currentTour
bool foundIn(Tour& currentTour, const int& i)
{
    return std::find(currentTour.cityVisitOrder.begin(), currentTour.cityVisitOrder.end(), i) != currentTour.cityVisitOrder.end();
}

// performs necessary actions to visit a city
void VisitCity(Tour& currentTour, std::set<int>& visitableCities, int i)
{
    currentTour.currentCity = i;
    currentTour.cityVisitOrder.push_back(i);

    visitableCities.erase(i);
}

// choose a random starting city
int ChooseStartingCity(int numcoordinates)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, numcoordinates - 1);

    int randomNum = dist(gen);
    return randomNum;
}

// implementation of Kruskals Algorithm to find the MST(Minimum Spanning tree
std::vector<Edge> findMST(const std::vector<Coordinate>& coordinates)
{
    int numCoordinates = coordinates.size();

    // create a vector of edges 
    std::vector<Edge> edges;

    // calculating the distances of all the edges and storing them in the "edges"
    for (int i = 0; i < numCoordinates; ++i)
    {
        for (int j = i + 1; j < numCoordinates; ++j)
        {
            double distance = calculateDistance(coordinates[i], coordinates[j]);
            Edge edge;
            edge.coord1 = coordinates[i];
            edge.coord2 = coordinates[j];
            edge.distance = distance;
            edges.push_back(edge);
        }
    }

    // sort the edges from shortest to longest
    std::sort(edges.begin(), edges.end(), [](const Edge& e1, const Edge& e2) 
    {
        return e1.distance < e2.distance;
    });

    // needed data structures
    std::vector<Edge> mst;
    
    std::vector<int> parent(numCoordinates);
    
    std::vector<int> rank(numCoordinates, 0);

    // initilization of "parent" vector
    for (int i = 0; i < numCoordinates; ++i)
        parent[i] = i;

    int numEdges = 0;
    int index = 0;

    // fill the mst by adding edges from shortest to longest
    while (numEdges < numCoordinates - 1)
    {
        Edge currentEdge = edges[index++];
        int root1 = findRoot(parent, currentEdge.coord1.id);
        int root2 = findRoot(parent, currentEdge.coord2.id);

        // check for a cycle
        if (root1 != root2)
        {
            // add the "currentEdge" to the "mst" vector of "Edge"s
            mst.push_back(currentEdge);
            unionNodes(parent, rank, root1, root2);
            numEdges++;
        }
    }

    return mst;
}

// find the odd degree vertices of the mst
std::vector<int> findOddDegreeVertices(const std::vector<Edge>& mst)
{
    std::vector<int> oddDegreeVertices;
    std::vector<int> degreeCount(mst.size() + 1, 0); // +1 to account for 1-based indexing of vertices

    for (const auto& edge : mst)
    {
        degreeCount[edge.coord1.id]++;
        degreeCount[edge.coord2.id]++;
    }

    for (int i = 1; i < degreeCount.size(); ++i)
    {
        if (degreeCount[i] % 2 != 0)
            oddDegreeVertices.push_back(i);
    }

    return oddDegreeVertices;
}

// find the shortest matching for the odd degree vertices to make them even degree
std::vector<Edge> findMinimumWeightPerfectMatching(const std::vector<Coordinate>& coordinates, std::vector<int>& oddDegreeVertices)
{    
    std::vector<Edge> matching;

    for (size_t i = 0; i < oddDegreeVertices.size(); ++i)
    {
        int u = oddDegreeVertices[i];
        double minDistance = std::numeric_limits<double>::max();
        int minIndex = -1;

        for (size_t j = i + 1; j < oddDegreeVertices.size(); ++j)
        {
            int v = oddDegreeVertices[j];
            double distance = calculateDistance(coordinates[u - 1], coordinates[v - 1]);

            if (distance < minDistance)
            {
                minDistance = distance;
                minIndex = j;
            }
        }

        if (minIndex != -1)
        {
            Edge edge;
            edge.coord1 = coordinates[u - 1];
            edge.coord2 = coordinates[oddDegreeVertices[minIndex] - 1];
            edge.distance = minDistance;
            matching.push_back(edge);
            oddDegreeVertices.erase(oddDegreeVertices.begin() + minIndex);
        }
    }

    return matching;
}

// combines the mst and best matching graphs back
std::vector<Edge> combineMSTAndMatching(const std::vector<Edge>& mst, const std::vector<Edge>& matching)
{
    std::vector<Edge> multigraph = mst;

    for (const auto& edge : matching)
    {
        multigraph.push_back(edge);
    }

    return multigraph;
}

// find the eulerian circuit 
std::vector<int> findEulerianCircuit(const std::vector<Edge>& multigraph)
{    
    std::vector<std::vector<int>> adjacencyList(multigraph.size() + 1);

    for (size_t i = 0; i < multigraph.size(); ++i)
    {
        int u = multigraph[i].coord1.id;
        int v = multigraph[i].coord2.id;

        adjacencyList[u].push_back(v);
        adjacencyList[v].push_back(u);
    }

    std::vector<int> eulerianCircuit;
    std::stack<int> stack;
    int currentVertex = multigraph[0].coord1.id;

    while (true)
    {
        if (!adjacencyList[currentVertex].empty())
        {
            stack.push(currentVertex);
            int nextVertex = adjacencyList[currentVertex].back();
            adjacencyList[currentVertex].pop_back();
            adjacencyList[nextVertex].erase(std::find(adjacencyList[nextVertex].begin(), adjacencyList[nextVertex].end(), currentVertex));
            currentVertex = nextVertex;
        }
        else if (!stack.empty())
        {
            eulerianCircuit.push_back(currentVertex);
            currentVertex = stack.top();
            stack.pop();
        }
        else
        {
            break;
        }
    }

    return eulerianCircuit;
}

// converts eulerian circuit into a hamiltonian circuit
std::vector<int> modifyEulerianCircuit(const std::vector<int>& eulerianCircuit)
{
    std::vector<int> hamiltonianCircuit;
    std::vector<bool> visited(eulerianCircuit.size() + 1, false); 

    for (int vertex : eulerianCircuit)
    {
        if (!visited[vertex])
        {
            hamiltonianCircuit.push_back(vertex);
            visited[vertex] = true;
        }
    }

    return hamiltonianCircuit;
}

// calculates the cost
double calculateHamiltonianCircuitCost(const std::vector<int>& hamiltonianCircuit, const std::vector<Coordinate>& coordinates)
{
    double totalCost = 0.0;

    for (size_t i = 0; i < hamiltonianCircuit.size(); ++i) {
        int currIndex = hamiltonianCircuit[i];
        int nextIndex = hamiltonianCircuit[(i + 1) % hamiltonianCircuit.size()];

        if (currIndex >= 0 && currIndex < coordinates.size() &&
            nextIndex >= 0 && nextIndex < coordinates.size()) {
            double distance = calculateDistance(coordinates[currIndex], coordinates[nextIndex]);
            totalCost += distance;
        }
        else {
            std::cout << "Out Of Range";
        }
    }

    return totalCost;
}

// swap two edges
void twoOptSwap(std::vector<int>& circuit, int i, int k)
{
    while (i < k) {
        std::swap(circuit[i], circuit[k]);
        i++;
        k--;
    }
}

// 2-opt optiomazition
void twoOpt(std::vector<int>& circuit, const std::vector<Coordinate>& coordinates)
{
    int numVertices = circuit.size();
    bool improvement = true;

    while (improvement) {
        improvement = false;

        for (int i = 0; i < numVertices - 2; i++) {
            for (int k = i + 2; k < numVertices; k++) {
                double distanceBefore = calculateDistance(coordinates[circuit[i]], coordinates[circuit[i + 1]])
                    + calculateDistance(coordinates[circuit[k]], coordinates[circuit[(k + 1) % numVertices]]);
                double distanceAfter = calculateDistance(coordinates[circuit[i]], coordinates[circuit[k]])
                    + calculateDistance(coordinates[circuit[i + 1]], coordinates[circuit[(k + 1) % numVertices]]);

                if (distanceAfter < distanceBefore) {
                    twoOptSwap(circuit, i + 1, k);
                    improvement = true;
                }
            }
        }
    }
}

// 3-opt optimization
void threeOpt(std::vector<int>& circuit, const std::vector<Coordinate>& coordinates)
{
    int numVertices = circuit.size();
    bool improvement = true;

    while (improvement) {
        improvement = false;

        for (int i = 0; i < numVertices - 4; i++) {
            for (int j = i + 2; j < numVertices - 2; j++) {
                for (int k = j + 2; k < numVertices; k++) {
                    // Case 1: 2-opt swap
                    double distanceBefore = calculateDistance(coordinates[circuit[i]], coordinates[circuit[i + 1]])
                        + calculateDistance(coordinates[circuit[j]], coordinates[circuit[j + 1]])
                        + calculateDistance(coordinates[circuit[k]], coordinates[circuit[(k + 1) % numVertices]]);
                    double distanceAfter = calculateDistance(coordinates[circuit[i]], coordinates[circuit[j]])
                        + calculateDistance(coordinates[circuit[i + 1]], coordinates[circuit[(j + 1) % numVertices]])
                        + calculateDistance(coordinates[circuit[k]], coordinates[circuit[(k + 1) % numVertices]]);

                    if (distanceAfter < distanceBefore) {
                        twoOptSwap(circuit, i + 1, j);
                        improvement = true;
                        continue;
                    }

                    // Case 2: 2-opt swap
                    distanceBefore = calculateDistance(coordinates[circuit[i]], coordinates[circuit[i + 1]])
                        + calculateDistance(coordinates[circuit[j]], coordinates[circuit[j + 1]])
                        + calculateDistance(coordinates[circuit[k]], coordinates[circuit[(k + 1) % numVertices]]);
                    distanceAfter = calculateDistance(coordinates[circuit[i]], coordinates[circuit[i + 1]])
                        + calculateDistance(coordinates[circuit[j]], coordinates[circuit[k]])
                        + calculateDistance(coordinates[circuit[j + 1]], coordinates[circuit[(k + 1) % numVertices]]);

                    if (distanceAfter < distanceBefore) {
                        twoOptSwap(circuit, j + 1, k);
                        improvement = true;
                        continue;
                    }

                    // Case 3: 2-opt swap
                    distanceBefore = calculateDistance(coordinates[circuit[i]], coordinates[circuit[i + 1]])
                        + calculateDistance(coordinates[circuit[j]], coordinates[circuit[j + 1]])
                        + calculateDistance(coordinates[circuit[k]], coordinates[circuit[(k + 1) % numVertices]]);
                    distanceAfter = calculateDistance(coordinates[circuit[i]], coordinates[circuit[i + 1]])
                        + calculateDistance(coordinates[circuit[j]], coordinates[circuit[k]])
                        + calculateDistance(coordinates[circuit[(k + 1) % numVertices]], coordinates[circuit[j + 1]]);

                    if (distanceAfter < distanceBefore) {
                        twoOptSwap(circuit, i + 1, k);
                        improvement = true;
                    }
                }
            }
        }
    }
}


int main() {
    // Choose Algoritm
    int userInput;    
    bool validInput = false;    

    //retrieved from the data source for uruguay case
    /*double optimalDistance = 79114;*/

    while (!validInput)
    {
        std::cout << "Do you wanna use Ant Colony Optimization algorithm or Christofides algorithm? Type 1 or 2 respectively: ";
        std::cin >> userInput;
        if (userInput == 1)
        { // ANT COLONY OPTIMIZATION ALGORITHM IS CHOSEN 
            validInput = true;
        }
        else if (userInput == 2)
        { // CHRISTOFIDES ALGORITHM IS CHOSEN
            validInput = true;
        }
        else
        {
            validInput = false;
            std::cout << "Wrong Number, please enter 1 or 2" << std::endl;
        }
    }

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

    double shortestDistance = 2 * std::pow(10, 8);

    if (userInput == 1)
    {
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


        for (int run = 0; run < 20000; run++)
        {
            int importanceEdgeRating = 0 + run/1000 * 2;
            int importanceInverseDistance = 15 - run/1000;
        
        

            // Tour is a struct that contains int(currnetCity), vector of ints(list of visited city ids in order in this case) and a double holding the total distance traveled in the tour
            Tour currentTour;

            // creating set for unvisited cities so that we can remove by value
            std::set<int> visitableCities;

            for (int i = 0; i < coordinates.size(); i++)
            {
                visitableCities.insert(i);
            }

            // choose random starting city
            int startCityIndex = ChooseStartingCity(numCoordinates);
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
                        double neighborRatingProb =  std::pow(edgeRatingMatrix[currentTour.currentCity][neighborProb.id] / totalRatingVisitableCities, importanceEdgeRating);
                        double neighborDistanceProb = std::pow(edgeInverseDistanceMatrix[currentTour.currentCity][neighborProb.id] / totalInverseDistanceVisitableCities, importanceInverseDistance);
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
                    // scale neighbor probability along with Total Probability
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

            if (run != 0 && run % 100 == 0)
            {
                std::cout << "hey";
            }
            // adjusting the edge rating matrix
            for (int i = 0; i < currentTour.cityVisitOrder.size() - 1; i++)
            {
                edgeRatingMatrix[currentTour.cityVisitOrder[i]][currentTour.cityVisitOrder[i + 1]] += std::pow(10, 2) / (currentTour.distanceTraveled/*-optimalDistance*/); //"-optimalDistance" can be used for more aggressive rewarding
                edgeRatingMatrix[currentTour.cityVisitOrder[i + 1]][currentTour.cityVisitOrder[i]] += std::pow(10, 2) / (currentTour.distanceTraveled/*-optimalDistance*/);
            }

            if (currentTour.distanceTraveled < shortestDistance)
            {
                shortestDistance = currentTour.distanceTraveled;
                bestTour = currentTour;
            }
            std::cout << "On Run: " << run << "\t";
            std::cout << "The distance was: " << currentTour.distanceTraveled << std::endl;
        }
        std::cout << "Shortest distance was: " << shortestDistance << std::endl;
        std::cout << "And the route was: " << std::endl;
        for (auto& cityId : bestTour.cityVisitOrder)
        {
            std::cout << cityId << "->";
        }
    }

    if (userInput == 2)
    {   // Christofides Algorithm

       // // Need to find MST(Minimum Spanning Tree)

        // Kruskals Algorithm to find the MST
        std::vector<Edge> mst = findMST(coordinates);

        // Find the odd degree vertices in "mst"
        std::vector<int> oddDegreeVertices = findOddDegreeVertices(mst);

        // Find the min-weight perfect matching for the odd degree vertices
        std::vector<Edge> matching = findMinimumWeightPerfectMatching(coordinates, oddDegreeVertices);

        // Combine the MST and the minimum-weight perfect matching together
        std::vector<Edge> multigraph = combineMSTAndMatching(mst, matching);

        // find the eulerian circuit
        std::vector<int> eulerianCircuit = findEulerianCircuit(multigraph);

        // Conver the eulerian circuit into a hamiltonian circuit
        std::vector<int> hamiltonianCircuit = modifyEulerianCircuit(eulerianCircuit);

        // Find the totalCost of the circuit
        double totalCost = calculateHamiltonianCircuitCost(hamiltonianCircuit, coordinates);

        // Further optimization of the circuit
        //while (calculateHamiltonianCircuitCost(hamiltonianCircuit, coordinates) > 89000)
        //{
            //twoOpt(hamiltonianCircuit, coordinates);
            //threeOpt(hamiltonianCircuit, coordinates);
        //}
        
        double optimizedCost = calculateHamiltonianCircuitCost(hamiltonianCircuit, coordinates);

        std::cout << "Total cost of the Hamiltonian circuit: " << optimizedCost << std::endl;

    }
    
    

    return 0;
}







