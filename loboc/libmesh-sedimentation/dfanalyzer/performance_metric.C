#include "performance_metric.h"

#include <ctime>

void PerformanceMetric::IdentifyStartTime() {
    time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];

  time (&rawtime);
  timeinfo = localtime(&rawtime);

  strftime(buffer,80,"%d.%m.%Y-%I:%M:%S",timeinfo);
  std::string str(buffer);

  this->startTime = str;
}

void PerformanceMetric::IdentifyEndTime() {
    time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];

  time (&rawtime);
  timeinfo = localtime(&rawtime);

  strftime(buffer,80,"%d.%m.%Y-%I:%M:%S",timeinfo);
  std::string str(buffer);

  this->endTime = str;
}