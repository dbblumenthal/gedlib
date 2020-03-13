#ifndef __IBD_DATASET_H__
#define __IBD_DATASET_H__

#include "Dataset.h"
#include "IBDGraph.h"

class IBDDataset : public Dataset<int, double, int>
{

protected:

  /**
   * Load a complete dataset, ie. a list of QAP Instances
   */
  void loadDS(const char* filename);

  /**
   * Load a QAP instance from a .dat file
   */
  void loadQAP(const char* filename);

public:

  IBDDataset(const char* filename);
  IBDDataset(): Dataset<int, double, int>(){}

};


#endif
