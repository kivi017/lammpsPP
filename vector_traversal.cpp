#include <iostream>
#include <vector>
#include <cmath>

int main(int argc, char const *argv[]) {
  std::vector<int> V;   //Creating a vector V of integer type

  for(int i=1; i<6; i++)
    V.push_back(3*i);   //Adding new element to the end of the vector V

  for(std::vector<int>::iterator it=V.begin(); it != V.end(); it++)   //Loop to iterate through the vector
    std::cout<<"\n"<<pow(*it-1,2)<<std::endl;

  return 0;
  }
