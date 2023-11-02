#include <iostream>
#include <stdlib.h>
#include <utility>
#include <cmath> // square function
#include <algorithm> // find function


int main () {
    int size;
    typedef std::pair<int,int> P;
    
    srand(time(nullptr));

    std::cout << "input graph size";
    std::cin >> size;

    P coords[size];
    float matrix[size][size] = {}; 
     

    for ( int i = 0, x, y; i < size; i++) { // setting coordinates for points
        y = rand() % 10000;
        x = rand() % 10000;
        P coordinates(x, y);
        std::cout << coordinates.first << " ";
        std::cout << coordinates.second << "\n";
        coords[i] = coordinates;
    }
    
    for ( int i = 0; i < size; i++){ // setting all [n][n] relations with -1, to not take them into account in next step
        matrix[i][i] = -1;
    }
    
    for ( int i = 0; i < size; i++) { // calculating distance between two points and puting it in matrix
        for ( int j = (i + 1), x; j < size; j++ ) {
            float distance = sqrt((coords[j].first - coords[i].first)*(coords[j].first - coords[i].first) + (coords[j].second - coords[i].second)*(coords[j].second - coords[i].second));
            matrix[i][j] = distance;
            matrix[j][i] = distance;
        }
    }
    
    
    
    int path[size], lenght = 0;
    path[0] = 0;
    
    for ( int i = 0, location = 0; i < size; i++){
        int current = 0;
        float currentValue = 100000000000000000;
        for ( int j = 0; j < size; j ++) { // finding closest point
            bool isPresent = std::find(path, path + size, j) != path + size;
            if (matrix[location][j] >= 0 && matrix[location][j] <= currentValue && !isPresent){
                current = j;
                currentValue = matrix[location][current];
                
            }
        }
        lenght = lenght + matrix[location][current];
        location = current;
        path[i+1] = location; // not at the start, to make sure last location is recorded
        
    }
    
    //showing results
    
    std::cout << "Lenght of path " << lenght << "\n Path:";
    
    for ( int i = 0; i < size; i++ ){
        std::cout << path[i] << " ";
    }
    return 0;
}