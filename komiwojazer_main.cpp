#include <iostream>
#include <stdlib.h>
#include <utility>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>

std::vector<std::pair<int, int>> GenerateRandomGraph(int size){
    typedef std::pair<int, int> P;
    std::vector<P> coords;
    srand(time(nullptr));
    
    for ( int i = 0, x, y; i < size; i++) { // setting coordinates for points
        y = rand() % 10000;
        x = rand() % 10000;
        P coordinates(x, y);
        std::cout << coordinates.first << " ";
        std::cout << coordinates.second << "\n";
        coords.push_back(coordinates);
    }
    return coords;
}
void SaveGraph(const char* filename, int size, std::vector<std::pair<int, int>> coords){
    std::ofstream outputFile(filename);
    if (outputFile.is_open()){
        outputFile << size << std::endl;
        for(int i = 0; i<size; i++){
            outputFile << i+1 << " " << coords[i].first << coords[i].second << std::endl;
        }

    } else {
        std::cout<<"Opening file failed"<<std::endl;
    }
}
std::vector<std::pair<int, int>> GetGraph(std::string filename) {
    std::ifstream inputFile(filename);
    std::vector<std::pair<int, int>> coords;
    int nr;
    std::string coord1, coord2;
    int size;
    if (inputFile.is_open()) {
        inputFile >> size;
        while(inputFile >> nr >> coord1 >> coord2){
            coords.push_back(std::pair<int,int>(std::stoi(coord1), std::stoi(coord2)));
        }
        return coords;
    } else {
        std::cout<<"Opening file failed"<<std::endl;
    } 
}
std::pair<std::vector<std::pair<int, int>>, std::string> NewInput(int argc, char *argv[]){
    std::string filename;
    
    typedef std::pair<int, int> P;
    std::vector<P> coords;
    int size;
    
    if(argc == 1){
        std::cout<<"Input graph size: "<<std::endl;
        std::cin>>size;
        coords = GenerateRandomGraph(size);  
        SaveGraph("graph.txt", size, coords);  
        filename = "graph.txt";
    } else if (argc == 2) {
        filename = argv[1];
        coords = GetGraph(filename);  
    }
    return std::pair<std::vector<std::pair<int, int>>, std::string> (coords, filename);
}

float GetDistance(std::vector<std::pair<int, int>> coords, int i, int j){
    return sqrt((coords[j].first - coords[i].first)*(coords[j].first - coords[i].first) + (coords[j].second - coords[i].second)*(coords[j].second - coords[i].second));
}
void SaveResult(int allDist, std::vector<int> path, std::string filename, int size){
    std::ofstream resultFile(filename);
    if(resultFile.is_open()) {
        resultFile << "Length of the path: "<< allDist << std::endl;
        resultFile << "Entire path: "<< std::endl;
        for ( int i = 0; i <= size; i++ ){
            std::cout << path[i] << " ";
        }
    } else {
        std::cout << "Error while opening the result file" <<std::endl;
    }

}
int main (int argc, char* argv[]) {
    std::string filename;
    std::pair<std::vector<std::pair<int, int>>, std::string> NewPair;
    std::vector<std::pair<int, int>> coords;
    NewPair = NewInput(argc, argv);
    coords = NewPair.first;
    filename = NewPair.second;
    int size = coords.size();
    typedef std::pair<int,int> P;
    std::vector<std::vector<float>> matrix;
    std::cout << size;
    for ( int i = 0; i < size; i++) { // calculating distance between two points and puting it in matrix
        matrix.push_back(std::vector<float>(size, -1));
        for (int j = 0; j < size; j++){
            if (i != j){
                float distance = GetDistance(coords, i, j);
            }
        }
    }
    for (int i =0; i < size; i++){
        std::cout << coords[i].first << " " << coords[i].second << std::endl;
    }
    std::cout<<"sa";
    std::vector<int> path;
    std::vector<int> visited(size, 0);
    int lenght = 0;
    path.push_back(0);
    visited[0] = 1;
    int dest;
    float allDist = 0.0;

    while(path.size() < size) {
        std::cout<<"sqdqw";
        int current = path.back();
        float currentValue = 100000000000000000;
        for (int j = 0; j < size; j ++) {
            if(current!=j && !visited[j] && matrix[current][j] >= 0 && matrix[current][j] <= currentValue){
                dest = j;
                currentValue = matrix[current][j];
            }
        }
        std::cout << "testCZZC";
        path.push_back(dest);
        std::cout << "teahjsgfkjashflksajfoisastCZZC";
        visited[dest] = 1;
        std::cout <<124124124;
        allDist += currentValue;
        std::cout << "test";

    }
    std::cout << "test`12312";
    path.push_back(0);
    allDist += matrix[0][dest];
    
    //showing results
    std::cout<<"fdre";
    std::cout << "Lenght of path " << allDist << "\nPath: ";
    
    for ( int i = 0; i <= size; i++ ){
        std::cout << path[i] << " ";
    }
    return 0;
}