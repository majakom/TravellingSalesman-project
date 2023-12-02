#include <iostream>
#include <stdlib.h>
#include <utility>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>


int main () {
    int size;
    typedef std::pair<int,int> P;
    
    srand(time(nullptr));
    std::cout << "File (0) input (1)";
    int input;
    std::cin >> input;
    std::vector<P> coords;
    std::vector<std::vector<float>> matrix;
    if (input){
        std::cout<<"Input graph size: "<<std::endl;
        std::cin >> size;

        for ( int i = 0, x, y; i < size; i++) { // setting coordinates for points
            y = rand() % 10000;
            x = rand() % 10000;
            P coordinates(x, y);
            std::cout << coordinates.first << " ";
            std::cout << coordinates.second << "\n";
            coords.push_back(coordinates);
        }
        std::ofstream outputFile("random.txt");
        if (outputFile.is_open()){
            outputFile << size << std::endl;
            for(int i = 0; i<size; i++){
                outputFile << i+1 << " " << coords[i].first <<" "<< coords[i].second << std::endl;
            }

        } else {
            std::cout<<"Opening file failed"<<std::endl;
        }
    }
    else{
        std::ifstream inputFile("test.txt");
        if(inputFile.is_open()){
            inputFile >> size;
            int edge_number, v1, v2;
            while (inputFile >> edge_number >> v1 >> v2){
                coords.push_back(std::pair<int, int>(v1, v2));   
            }  
        }else{
            std::cout << "Error. Couldn't open file: " << "test.txt" << std::endl;
        }
    }
    
    
    for ( int i = 0; i < size; i++) { // calculating distance between two points and puting it in matrix
        matrix.push_back(std::vector<float>(size, -1));
        for (int j = 0; j < size; j++){
            if (i != j){
                float distance = sqrt((coords[j].first - coords[i].first)*(coords[j].first - coords[i].first) + (coords[j].second - coords[i].second)*(coords[j].second - coords[i].second));
                matrix[i][j] = distance;
            }
        }
    }
    
    
    std::vector<int> path;
    std::vector<int> visited(size, 0);
    int lenght = 0;
    path.push_back(0);
    visited[0] = 1;
    int dest;
    float allDist = 0.0;

    while(path.size() < size) {
        int current = path.back();
        float currentValue = 100000000000000000;
        for (int j = 0; j < size; j ++) {
            if(current!=j && !visited[j] && matrix[current][j] >= 0 && matrix[current][j] <= currentValue){
                dest = j;
                currentValue = matrix[current][j];
            }
        }
        path.push_back(dest);
        visited[dest] = 1;
        allDist += currentValue;

    }
    path.push_back(0);
    allDist += matrix[0][dest];
    
    //showing results
    std::ofstream outputFile("out.txt");
    if (outputFile.is_open()){
        outputFile << "Length of path: " << allDist << std::endl;
        for(int i = 0; i<=size; i++){
            outputFile << i+1 << " " << path[i]+1 << std::endl;
        }

    } else {
        std::cout<<"Opening file failed"<<std::endl;
    }
    
    std::cout << "Lenght of path " << allDist << "\nPath: ";
    
    for ( int i = 0; i <= size; i++ ){
        std::cout << path[i] + 1 << " ";
    }
    return 0;
}