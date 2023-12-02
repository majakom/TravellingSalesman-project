#include <iostream>
#include <stdlib.h>
#include <utility>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>

std::vector<std::pair<int, int>> ChooseInput(int argc, char* argv[]);
std::vector<std::pair<int, int>> GenerateRandomGraph(int size);
std::vector<std::pair<int, int>> GetGraphFromFile(std::string file);
std::vector<std::vector<double>> Matrix(std::vector<std::pair<int, int>> coords, int size);
std::pair<double, std::vector<int>> greedyAlgorithmTSP(int size, std::vector<std::vector<double>> matrix);
void OutputToTerminal(double allDistance, std::vector<int> path, int size);
void SaveOutput(double allDistance, std::vector<int> path);
double GetDistance(std::vector<std::pair<int, int>> coords, int i, int j);

int main(int argc, char *argv[]){
    std::vector<std::pair<int, int>> coords;
    std::vector<std::vector<double>> matrix;
    std::pair<double, std::vector<int>> output;
    std::vector<int> path;
    double allDistance;
    
    
    coords = ChooseInput(argc, argv);
    int size = coords.size();
    matrix = Matrix(coords, size);
    output = greedyAlgorithmTSP(size, matrix);

    SaveOutput(output.first, output.second);
    OutputToTerminal(output.first, output.second, size);
    return 0;
}

std::vector<std::pair<int, int>> ChooseInput(int argc, char* argv[]){
    std::string file;
    std::vector<std::pair<int, int>> coords;
    int size;

    if (argc==2){
        file = argv[1];
        coords = GetGraphFromFile(file);
    } else if (argc == 1) {
        std::cout<<"Input graph size: "<<std::endl;
        std::cin>>size;
        coords = GenerateRandomGraph(size);
    }
    return coords;
}
std::vector<std::pair<int, int>> GenerateRandomGraph(int size) {
    std::ofstream outputFile("random.txt");
    std::vector<std::pair<int, int>> coords;
    std::pair<int, int> newPair;
    int coord1, coord2;

    srand(time(nullptr));
    

    for (int i = 0; i < size; i++) { // setting coordinates for points
        coord1 = rand() % 10000;
        coord2 = rand() % 10000;
        newPair.first = coord1;
        newPair.second = coord2;
        std::cout << newPair.first << " ";
        std::cout << newPair.second << "\n";
        coords.push_back(newPair);
    }
    if (outputFile.is_open()) {
        outputFile << size << std::endl;
        for(int i = 0; i<size; i++){
            outputFile << i+1 << " " << coords[i].first <<" "<< coords[i].second << std::endl;
        }
    } else {
        std::cout<<"Opening file for a random graph input failed"<<std::endl;
    }
    return coords;
}
std::vector<std::pair<int, int>> GetGraphFromFile(std::string file) {
    std::ifstream inputFile(file);
    std::vector<std::pair<int, int>> coords;
    std::string coord1, coord2;
    int size, nr; 
    if (inputFile.is_open()) {
        inputFile >> size;
        while(inputFile >> nr >> coord1 >> coord2){
            coords.push_back(std::pair<int,int>(std::stoi(coord1), std::stoi(coord2)));
        }
        return coords;
    } else {
        std::cout<<"Opening file with coordinates failed"<<std::endl;
    }
    return coords; 

}
double GetDistance(std::vector<std::pair<int, int>> coords, int i, int j){
    return sqrt((coords[j].first - coords[i].first)*(coords[j].first - coords[i].first) + (coords[j].second - coords[i].second)*(coords[j].second - coords[i].second));
}
std::vector<std::vector<double>> Matrix(std::vector<std::pair<int, int>> coords, int size){
    std::vector<std::vector<double>> matrix;
    for ( int i = 0; i < size; i++) {
        matrix.push_back(std::vector<double>(size, -1));
        for (int j = 0; j < size; j++){
            if (i != j){
                matrix[i][j] = GetDistance(coords, i, j);
            }
        }
    }
    return matrix;
}
std::pair<double, std::vector<int>> greedyAlgorithmTSP(int size, std::vector<std::vector<double>> matrix){
    std::vector<int> path;
    std::vector<int> visited(size, 0);
    double allDistance = 0.0;
    int currentVertex, destination;
    double currentValue;

    visited[0] = 1;
    path.push_back(0);
    
    while(path.size()<size){
        currentValue = 100000000000000000;
        currentVertex = path.back();
        for (int i = 0; i < size; i++){
            if(currentVertex != i && !visited[i] && matrix[currentVertex][i] < currentValue) {
                currentValue = matrix[currentVertex][i];
                destination = i;
            }
        }
        path.push_back(destination);
        visited[destination] = 1;
        allDistance += currentValue;
    }
    path.push_back(0);
    allDistance += matrix[0][destination];

    return std::pair<double, std::vector<int>>(allDistance, path);
}
void SaveOutput(double allDistance, std::vector<int> path) {
    std::ofstream resultFile("output.txt");
    if (resultFile.is_open()) {
        resultFile << "Distance: "<< allDistance <<std::endl;
        for (int i = 0; i < path.size(); i++){
            resultFile << path[i] + 1 << "\n";
        }
        resultFile << std::endl;
    } else {
        std::cout<<"Could not open result file"<<std::endl;
    }
}
void OutputToTerminal(double allDistance, std::vector<int> path, int size) {
    std::cout << "Length of path: " << allDistance << std::endl << "Path: ";
    for ( int i = 0; i <= size; i++ ){
        std::cout << path[i] + 1 << " ";
    }
}