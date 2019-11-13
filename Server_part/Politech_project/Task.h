#pragma once
#include <map>
#include <string>
#include "Point.h"
class Task
{
public:
    Task(){};
    ~Task(){};
    virtual void nextStep()=0;
    virtual double getTime()=0;
    virtual std::map<std::string, Point> getPoints()=0;
    virtual void parseParam(std::map<std::string, std::string> mapParam)=0;
};

