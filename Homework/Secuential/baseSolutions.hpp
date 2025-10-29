//
//  baseSolutions.hpp
//  KnightTour
//
//  Created by Serguei Diaz on 22.12.2024.
//

#ifndef BASESOLUTIONS_HPP
#define BASESOLUTIONS_HPP

#include <iostream>
#include <map>
#include <vector>
#include <utility>

using namespace std;
using SolutionMatrix = vector<std::vector<int>>;

//map<pair<int, int>, SolutionMatrix> knightTourSolutions;

void initBaseSolutionsDoubleLoop();

void initBaseSolutionsOpen();

void initBaseSolutionsStretchedV();

void initBaseSolutionsStretchedH();

SolutionMatrix getStretchedSolutionV(int n, int m);

SolutionMatrix getStretchedSolutionH(int n, int m);

SolutionMatrix getDoubleLoopSolution(int n, int m);

SolutionMatrix getOpenSolution(int n, int m);

#endif

