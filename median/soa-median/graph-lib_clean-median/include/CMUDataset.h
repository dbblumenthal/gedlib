#ifndef __CMU_DATASET_H__
#define __CMU_DATASET_H__

#include "Dataset.h"
#include "CMUGraph.h"

class CMUDataset : public Dataset<CMUPoint, double, int>
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

  CMUDataset(const char* filename);

};


#endif //__CMU_DATASET_H__
